package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.cmdLine.Argument;

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
    @Argument(fullName="target_titv", shortName="titv", doc="The target Ti/Tv ratio towards which to optimize. (~~2.2 for whole genome experiments)", required=true)
    private double TARGET_TITV = 2.12;
    //@Argument(fullName="filter_output", shortName="filter", doc="If specified the optimizer will not only update the QUAL field of the output VCF file but will also filter the variants", required=false)
    //private boolean FILTER_OUTPUT = false;
    @Argument(fullName="ignore_input_filters", shortName="ignoreFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_INPUT_FILTERS = false;
    @Argument(fullName="exclude_annotation", shortName="exclude", doc="The names of the annotations which should be excluded from the calculations", required=false)
    private String[] EXCLUDED_ANNOTATIONS = null;
    @Argument(fullName="force_annotation", shortName="force", doc="The names of the annotations which should be forced into the calculations even if they aren't present in every variant", required=false)
    private String[] FORCED_ANNOTATIONS = null;
    @Argument(fullName="output", shortName="output", doc="The output file name", required=false)
    private String OUTPUT_FILE = "optimizer.data";
    @Argument(fullName="numGaussians", shortName="nG", doc="The number of Gaussians to be used in the Gaussian mixture model", required=false)
    private int NUM_GAUSSIANS = 32;
    @Argument(fullName="numIterations", shortName="nI", doc="The number of iterations to be performed in the Gaussian mixture model", required=false)
    private int NUM_ITERATIONS = 5; //BUGBUG: should automatically decided when to stop by looking at how entropy changes with each iteration


    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final ExpandingArrayList<String> annotationKeys = new ExpandingArrayList<String>();
    private boolean firstVariant = true;
    private int numAnnotations = 0;
    private static final double INFINITE_ANNOTATION_VALUE = 10000.0;

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

            if( vc != null && vc.isSNP() && (IGNORE_INPUT_FILTERS || !vc.isFiltered()) ) { 
                if( firstVariant ) { // This is the first variant encountered so set up the list of annotations
                    annotationKeys.addAll( vc.getAttributes().keySet() );
                    if( annotationKeys.contains("ID") ) { annotationKeys.remove("ID"); } // ID field is added to the vc's INFO field?
                    if( annotationKeys.contains("DB") ) { annotationKeys.remove("DB"); }
                    if( EXCLUDED_ANNOTATIONS != null ) {
                        for( final String excludedAnnotation : EXCLUDED_ANNOTATIONS ) {
                            if( annotationKeys.contains( excludedAnnotation ) ) { annotationKeys.remove( excludedAnnotation ); }
                        }
                    }
                    if( FORCED_ANNOTATIONS != null ) {
                        for( final String forcedAnnotation : FORCED_ANNOTATIONS ) {
                            if( !annotationKeys.contains( forcedAnnotation ) ) { annotationKeys.add( forcedAnnotation ); }
                        }
                    }
                    numAnnotations = annotationKeys.size() + 1; // +1 for variant quality ("QUAL")
                    annotationValues = new double[numAnnotations];
                    firstVariant = false;
                }

                int iii = 0;
                for( final String key : annotationKeys ) {

                    double value = 0.0;
                    try {
                        value = Double.parseDouble( (String)vc.getAttribute( key, "0.0" ) );
                        if( Double.isInfinite(value) ) {
                            value = ( value > 0 ? 1.0 : -1.0 ) * INFINITE_ANNOTATION_VALUE;
                        }
                    } catch( NumberFormatException e ) {
                        // do nothing, default value is 0.0,
                    }
                    annotationValues[iii++] = value;
                }

                // Variant quality ("QUAL") is not in the list of annotations, but is useful so add it here.
                annotationValues[iii] = vc.getPhredScaledQual();

                VariantDatum variantDatum = new VariantDatum();
                variantDatum.annotations = annotationValues;
                variantDatum.isTransition = vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                variantDatum.isKnown = !vc.getAttribute("ID").equals(".");
                variantDatum.isFiltered = vc.isFiltered(); // BUGBUG: This field won't be needed in the final version.
                mapList.add( variantDatum );
            }
        }

        return mapList;
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

        final VariantDataManager dataManager = new VariantDataManager( reduceSum, annotationKeys);
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        logger.info( "There are " + dataManager.numVariants + " variants and " + dataManager.numAnnotations + " annotations.");
        logger.info( "The annotations are: " + annotationKeys + " and QUAL." );

        dataManager.normalizeData(); // Each data point is now [ (x - mean) / standard deviation ]
        
        // Create either the Gaussian Mixture Model or the Nearest Neighbors model and run it
        final VariantOptimizationModel gmm = new VariantGaussianMixtureModel( dataManager, TARGET_TITV, NUM_GAUSSIANS, NUM_ITERATIONS  );
        gmm.run( OUTPUT_FILE );
    }

}
