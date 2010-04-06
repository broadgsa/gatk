package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
 * Applies calibrated variant cluster parameters to variant calls to produce an accurate and informative variant quality score
 *
 * @author rpoplin
 * @since Mar 17, 2010
 *
 * @help.summary Applies calibrated variant cluster parameters to variant calls to produce an accurate and informative variant quality score
 */

public class ApplyVariantClustersWalker extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="target_titv", shortName="titv", doc="The target Ti/Tv ratio towards which to optimize. (~~2.1 for whole genome experiments)", required=true)
    private double TARGET_TITV = 2.1;
    @Argument(fullName="backOff", shortName="backOff", doc="The Gaussian back off factor, used to prevent overfitting by spreading out the Gaussians.", required=false)
    private double BACKOFF_FACTOR = 1.0;
    @Argument(fullName="desired_num_variants", shortName="dV", doc="The desired number of variants to keep in a theoretically filtered set", required=false)
    private int DESIRED_NUM_VARIANTS = 0;
    @Argument(fullName="ignore_all_input_filters", shortName="ignoreAllFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_ALL_INPUT_FILTERS = false;
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the optimizer will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="known_prior", shortName="knownPrior", doc="A prior on the quality of known variants, a phred scaled probability of being true. Setting to 0.0 means unused.", required=false)
    private double KNOWN_VAR_QUAL_PRIOR = 0.0;
    @Argument(fullName="output_prefix", shortName="output", doc="The prefix added to output VCF file name and optimization curve pdf file name", required=false)
    private String OUTPUT_PREFIX = "optimizer";
    @Argument(fullName="clusterFile", shortName="clusterFile", doc="The output cluster file", required=true)
    private String CLUSTER_FILENAME = "optimizer.cluster";
    //@Argument(fullName = "optimization_model", shortName = "om", doc = "Optimization calculation model to employ -- GAUSSIAN_MIXTURE_MODEL is currently the default, while K_NEAREST_NEIGHBORS is also available for small callsets.", required = false)
    private VariantOptimizationModel.Model OPTIMIZATION_MODEL = VariantOptimizationModel.Model.GAUSSIAN_MIXTURE_MODEL;
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantGaussianMixtureModel theModel = null;
    private VCFWriter vcfWriter;
    private Set<String> ignoreInputFilterSet = null;


    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet = new TreeSet<String>(Arrays.asList(IGNORE_INPUT_FILTERS));
        }

        switch (OPTIMIZATION_MODEL) {
            case GAUSSIAN_MIXTURE_MODEL:
                theModel = new VariantGaussianMixtureModel( TARGET_TITV, CLUSTER_FILENAME, BACKOFF_FACTOR );
                break;
            //case K_NEAREST_NEIGHBORS:
            //    theModel = new VariantNearestNeighborsModel( dataManager, TARGET_TITV, NUM_KNN );
            //    break;
            default:
                throw new StingException( "Variant Optimization Model is unrecognized. Implemented options are GAUSSIAN_MIXTURE_MODEL and K_NEAREST_NEIGHBORS" );
        }


        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("OQ", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "The original variant quality score"));
        hInfo.add(new VCFHeaderLine("source", "VariantOptimizer"));
        vcfWriter = new VCFWriter( new File(OUTPUT_PREFIX + ".vcf") );
        final TreeSet<String> samples = new TreeSet<String>();
        SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, new HashMap<Pair<String, String>, String>());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
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


        for( final GATKFeature feature : tracker.getAllRods() ) {
            Object rod = feature.getUnderlyingObject();
            if( rod != null && rod instanceof RodVCF ) {
                final RodVCF rodVCF = ((RodVCF) rod);
                //BUGBUG: figure out how to make this use VariantContext to be consistent with other VariantOptimizer walkers
                //          need to convert vc into VCFRecord to write it out?
                if( rodVCF.isSNP() &&
                        (!rodVCF.isFiltered() || IGNORE_ALL_INPUT_FILTERS || (ignoreInputFilterSet != null && ignoreInputFilterSet.containsAll(Arrays.asList(rodVCF.getFilteringCodes())))) ) {

                    final VariantDatum variantDatum = new VariantDatum();
                    variantDatum.isTransition = BaseUtils.isTransition((byte)rodVCF.getAlternativeBaseForSNP(), (byte)rodVCF.getReferenceForSNP()); //vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                    variantDatum.isKnown = !rodVCF.isNovel(); //!vc.getAttribute("ID").equals(".");
                    int numHet = 0;
                    int numHom = 0;
                    for( final VCFGenotypeRecord rec : rodVCF.getVCFGenotypeRecords() ) {
                        if( rec.isHet() ) { numHet++; }
                        else if( rec.isHom() ) { numHom++; }
                    }
                    variantDatum.isHet = numHet > numHom; //vc.getHetCount() > vc.getHomVarCount(); // BUGBUG: what to do here for multi sample calls?

                    final double pTrue = theModel.evaluateVariant( rodVCF.getInfoValues(), rodVCF.getQual(), variantDatum.isHet );
                    final double recalQual = QualityUtils.phredScaleErrorRate( Math.max( 1.0 - pTrue, 0.000000001) );                             

                    if( variantDatum.isKnown && KNOWN_VAR_QUAL_PRIOR > 0.1 ) {
                        variantDatum.qual = 0.5 * recalQual + 0.5 * KNOWN_VAR_QUAL_PRIOR;
                    } else {
                        variantDatum.qual = recalQual;
                    }
                    mapList.add( variantDatum );

                    rodVCF.mCurrentRecord.addInfoField("OQ", ((Double)rodVCF.getQual()).toString() );
                    rodVCF.mCurrentRecord.setQual( variantDatum.qual );
                    rodVCF.mCurrentRecord.setFilterString(VCFRecord.UNFILTERED);
                    vcfWriter.addRecord( rodVCF.mCurrentRecord );


                } else { // not a SNP or is filtered so just dump it out to the VCF file
                    vcfWriter.addRecord( rodVCF.mCurrentRecord );
                }
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

        vcfWriter.close();
        
        final VariantDataManager dataManager = new VariantDataManager( reduceSum, theModel.dataManager.annotationKeys );
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        theModel.outputOptimizationCurve( dataManager.data, OUTPUT_PREFIX, DESIRED_NUM_VARIANTS );

        // Execute Rscript command to plot the optimization curve
        // Print out the command line to make it clear to the user what is being executed and how one might modify it
        final String rScriptCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_OptimizationCurve.R" + " " + OUTPUT_PREFIX + ".dat" + " " + TARGET_TITV;
        System.out.println( rScriptCommandLine );

        // Execute the RScript command to plot the table of truth values
        try {
            Runtime.getRuntime().exec( rScriptCommandLine );
        } catch ( IOException e ) {
            throw new StingException( "Unable to execute RScript command: " + rScriptCommandLine );
        }
    }
}

