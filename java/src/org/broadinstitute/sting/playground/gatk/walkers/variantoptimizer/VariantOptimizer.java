package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.PrintStream;

/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing a callset that is optimized for
 *  an expected transition / transversion ratio.
 *
 * @author rpoplin
 * @since Feb 11, 2010
 *
 * @help.summary Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing an optimized callset
 */

public class VariantOptimizer extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName = "target_titv", shortName="titv", doc="The target Ti/Tv ratio to optimize towards. (~~2.2 for whole genome experiments)", required=true)
    public double TARGET_TITV = 2.2;
    @Argument(fullName = "filter_output", shortName="filter", doc="If specified the optimizer will not only update the QUAL field of the output VCF file but will also filter the variants", required=false)
    public boolean FILTER_OUTPUT = false;
    
    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final ExpandingArrayList<String> annotationKeys = new ExpandingArrayList<String>();
    private boolean firstVariant = true;
    private int numAnnotations = 0;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public ExpandingArrayList<VariantDatum> map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final ExpandingArrayList<VariantDatum> mapList = new ExpandingArrayList<VariantDatum>();

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return mapList;
        }

        double annotationValues[] = new double[numAnnotations];

        for( final VariantContext vc : tracker.getAllVariantContexts(null, context.getLocation(), false, false) ) {
            if( vc != null && vc.isSNP() ) {
                if( firstVariant ) { // This is the first variant encountered so set up the list of annotations
                    annotationKeys.addAll( vc.getAttributes().keySet() );
                    if( annotationKeys.contains("ID") ) { annotationKeys.remove("ID"); } // ID field is added to the vc's INFO field??
                    if( annotationKeys.contains("DB") ) { annotationKeys.remove("DB"); }
                    if( annotationKeys.contains("Dels") ) { annotationKeys.remove("Dels"); }
                    if( annotationKeys.contains("AN") ) { annotationKeys.remove("AN"); }
                    numAnnotations = annotationKeys.size() + 1; // +1 for variant quality ("QUAL")
                    annotationValues = new double[numAnnotations];
                    firstVariant = false;
                }

                //BUGBUG: for now only using the novel SNPs
                if( vc.getAttribute("ID").equals(".") ) {
                    int iii = 0;
                    for( final String key : annotationKeys ) {

                        double value = 0.0f;
                        try {
                            value = Double.parseDouble( (String)vc.getAttribute( key, "0.0" ) );
                        } catch( NumberFormatException e ) {
                            // do nothing, default value is 0.0f, annotations with zero variance will be ignored later
                        }
                        annotationValues[iii++] = value;
                    }

                    // Variant quality ("QUAL") is not in the list of annotations
                    annotationValues[iii] = vc.getPhredScaledQual();

                    VariantDatum variantDatum = new VariantDatum();
                    variantDatum.annotations = annotationValues;
                    variantDatum.isTransition = vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                    variantDatum.isKnown = !vc.getAttribute("ID").equals(".");
                    variantDatum.isFiltered = vc.isFiltered();
                    mapList.add( variantDatum );
                }
            }
        }

        return mapList; // This value isn't actually used for anything
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public ExpandingArrayList<VariantDatum> reduceInit() {
        return new ExpandingArrayList<VariantDatum>();
    }

    public ExpandingArrayList<VariantDatum> reduce( final ExpandingArrayList<VariantDatum> mapValue, final ExpandingArrayList<VariantDatum> reduceSum ) {
        reduceSum.addAll( mapValue );
        return reduceSum;
    }

    public void onTraversalDone( ExpandingArrayList<VariantDatum> reduceSum ) {

        final VariantDataManager dataManager = new VariantDataManager( reduceSum );
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        logger.info( "There are " + dataManager.numVariants + " variants and " + dataManager.numAnnotations + " annotations.");
        logger.info( "The annotations are: " + annotationKeys + " and QUAL." );

        dataManager.normalizeData(); // Each data point is now [ (x - mean) / standard deviation ]

        final VariantOptimizationModel gmm = new VariantGaussianMixtureModel( dataManager, TARGET_TITV );
        final double[][] p = gmm.run();
        final int numIterations = 16;

        try {
            final PrintStream out = new PrintStream("gmm512x16norm.data");
            for(int iii = 0; iii < dataManager.numVariants; iii++) {
                for( int ttt = 0; ttt < numIterations; ttt++ ) {
                    out.print(p[ttt][iii] + "\t");
                }
                out.println((dataManager.data[iii].isTransition ? 1 : 0)
                        + "\t" + (dataManager.data[iii].isKnown? 1 : 0)
                        + "\t" + (dataManager.data[iii].isFiltered ? 1 : 0)
                );
            }


        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

}
