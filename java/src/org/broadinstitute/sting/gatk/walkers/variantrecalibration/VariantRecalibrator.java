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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriterImpl;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Applies calibrated variant cluster parameters to variant calls to produce an accurate and informative variant quality score
 *
 * @author rpoplin
 * @since Mar 17, 2010
 *
 * @help.summary Applies calibrated variant cluster parameters to variant calls to produce an accurate and informative variant quality score
 */

public class VariantRecalibrator extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="target_titv", shortName="titv", doc="The expected Ti/Tv ratio to display on optimization curve output figures. (~~2.1 for whole genome experiments)", required=false)
    private double TARGET_TITV = 2.1;
    @Argument(fullName="backOff", shortName="backOff", doc="The Gaussian back off factor, used to prevent overfitting by spreading out the Gaussians.", required=false)
    private double BACKOFF_FACTOR = 1.3;
    @Argument(fullName="desired_num_variants", shortName="dV", doc="The desired number of variants to keep in a theoretically filtered set", required=false)
    private int DESIRED_NUM_VARIANTS = 0;
    @Argument(fullName="ignore_all_input_filters", shortName="ignoreAllFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_ALL_INPUT_FILTERS = false;
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the optimizer will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="hapmap_prior", shortName="hapmapPrior", doc="A prior on the quality of known variants, a phred scaled probability of being true.", required=false)
    private double HAPMAP_QUAL_PRIOR = 15.0;
    @Argument(fullName="known_prior", shortName="knownPrior", doc="A prior on the quality of known variants, a phred scaled probability of being true.", required=false)
    private double KNOWN_QUAL_PRIOR = 10.0;
    @Argument(fullName="novel_prior", shortName="novelPrior", doc="A prior on the quality of novel variants, a phred scaled probability of being true.", required=false)
    private double NOVEL_QUAL_PRIOR = 2.0;
    @Argument(fullName="output_prefix", shortName="output", doc="The prefix added to output VCF file name and optimization curve pdf file name", required=false)
    private String OUTPUT_PREFIX = "optimizer";
    @Argument(fullName="clusterFile", shortName="clusterFile", doc="The output cluster file", required=true)
    private String CLUSTER_FILENAME = "optimizer.cluster";
    @Argument(fullName="FDRtranche", shortName="tranche", doc="The levels of novel false discovery rate (FDR, implied by ti/tv) at which to slice the data. (in percent, that is 1.0 for 1 percent)", required=false)
    private Double[] FDR_TRANCHES = null;
    //@Argument(fullName = "optimization_model", shortName = "om", doc = "Optimization calculation model to employ -- GAUSSIAN_MIXTURE_MODEL is currently the default, while K_NEAREST_NEIGHBORS is also available for small callsets.", required = false)
    private VariantOptimizationModel.Model OPTIMIZATION_MODEL = VariantOptimizationModel.Model.GAUSSIAN_MIXTURE_MODEL;
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is maybe /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";
    @Argument(fullName="quality_step", shortName="qStep", doc="Resolution in QUAL units for optimization and tranche calculations", required=false)
    private double QUAL_STEP = 0.1;
    @Argument(fullName="singleton_fp_rate", shortName="fp_rate", doc="Prior expectation that a singleton call would be a FP", required=false)
    private double SINGLETON_FP_RATE = 0.5;
    @Argument(fullName="quality_scale_factor", shortName="qScale", doc="Multiply all final quality scores by this value. Needed to normalize the quality scores.", required=false)
    private double QUALITY_SCALE_FACTOR = 20.0;

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
                if ( SINGLETON_FP_RATE != -1 )
                    theModel.setSingletonFPRate(SINGLETON_FP_RATE);
                break;
            //case K_NEAREST_NEIGHBORS:
            //    theModel = new VariantNearestNeighborsModel( dataManager, TARGET_TITV, NUM_KNN );
            //    break;
            default:
                throw new StingException( "Variant Optimization Model is unrecognized. Implemented options are GAUSSIAN_MIXTURE_MODEL and K_NEAREST_NEIGHBORS" );
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        final TreeSet<String> samples = new TreeSet<String>();
        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("OQ", 1, VCFHeaderLineType.Float, "The original variant quality score"));
        hInfo.add(new VCFHeaderLine("source", "VariantOptimizer"));
        samples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit()));

        vcfWriter = new VCFWriterImpl( new File(OUTPUT_PREFIX + ".vcf") );
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);

        boolean foundDBSNP = false;
        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if( rod.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                foundDBSNP = true;
            }
        }

        if(!foundDBSNP) {
            throw new StingException("dbSNP track is required. This calculation is critically dependent on being able to distinguish known and novel sites.");
        }

        // Set up default values for the FDR tranches if necessary
        if( FDR_TRANCHES == null ) {
            FDR_TRANCHES = new Double[6];
            FDR_TRANCHES[0] = 0.1;
            FDR_TRANCHES[1] = 1.0;
            FDR_TRANCHES[2] = 5.0;
            FDR_TRANCHES[3] = 12.5;
            FDR_TRANCHES[4] = 25.0;
            FDR_TRANCHES[5] = 35.0;
        }
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

        for( final VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {
            if( vc != null && !vc.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) && vc.isSNP() ) {
                if( !vc.isFiltered() || IGNORE_ALL_INPUT_FILTERS || (ignoreInputFilterSet != null && ignoreInputFilterSet.containsAll(vc.getFilters())) ) {
                    final VariantDatum variantDatum = new VariantDatum();
                    variantDatum.isTransition = VariantContextUtils.getSNPSubstitutionType(vc).compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;

                    final DbSNPFeature dbsnp = DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
                    variantDatum.isKnown = dbsnp != null;
                    variantDatum.alleleCount = vc.getChromosomeCount(vc.getAlternateAllele(0)); // BUGBUG: assumes file has genotypes. Also, what to do about tri-allelic sites?

                    final double acPrior = theModel.getAlleleCountPrior( variantDatum.alleleCount );
                    double knownPrior = variantDatum.isKnown ? QualityUtils.qualToProb(KNOWN_QUAL_PRIOR) : QualityUtils.qualToProb(NOVEL_QUAL_PRIOR);
                    if(variantDatum.isKnown && DbSNPHelper.isHapmap( dbsnp ) ) {
                        knownPrior = QualityUtils.qualToProb(HAPMAP_QUAL_PRIOR);
                    }
                    final double totalPrior = 1.0 - ((1.0 - acPrior) * (1.0 - knownPrior));

                    if( MathUtils.compareDoubles(totalPrior, 1.0, 1E-8) == 0 || MathUtils.compareDoubles(totalPrior, 0.0, 1E-8) == 0 ) {
                        throw new StingException("Something is wrong with the prior that was entered by the user:  Prior = " + totalPrior); // TODO - fix this up later
                    }

                    final double pVar = theModel.evaluateVariant( vc );

                    final double lod = (Math.log10(totalPrior) + Math.log10(pVar)) - ((Math.log10(1.0 - totalPrior)) + Math.log10(1.0));
                    variantDatum.qual = Math.abs( QUALITY_SCALE_FACTOR * QualityUtils.lodToPhredScaleErrorRate(lod) );

                    mapList.add( variantDatum );

                    Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
                    attrs.put("OQ", String.format("%.2f", ((Double)vc.getPhredScaledQual())));
                    Set<String> filters = new HashSet<String>();
                    filters.add(VCFConstants.PASSES_FILTERS_v4);
                    VariantContext newVC = new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.getGenotypes(), variantDatum.qual / 10.0, filters, attrs);

                    vcfWriter.add( newVC, ref.getBase() );

                } else { // not a SNP or is filtered so just dump it out to the VCF file
                    vcfWriter.add( vc, ref.getBase() ); 
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

        theModel.outputOptimizationCurve( dataManager.data, OUTPUT_PREFIX, DESIRED_NUM_VARIANTS, FDR_TRANCHES, QUAL_STEP );

        // Execute Rscript command to plot the optimization curve
        // Print out the command line to make it clear to the user what is being executed and how one might modify it
        final String rScriptOptimizationCurveCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_OptimizationCurve.R" + " " + OUTPUT_PREFIX + ".dat" + " " + TARGET_TITV;
        final String rScriptTranchesCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Tranches.R" + " " + OUTPUT_PREFIX + ".dat" + " " + TARGET_TITV;
        logger.info( rScriptOptimizationCurveCommandLine );
        logger.info( rScriptTranchesCommandLine );

        // Execute the RScript command to plot the table of truth values
        try {
            Runtime.getRuntime().exec( rScriptOptimizationCurveCommandLine );
            Runtime.getRuntime().exec( rScriptTranchesCommandLine );
        } catch ( IOException e ) {
            Utils.warnUser("Unable to execute the RScript command.  While not critical to the calculations themselves, the script outputs a report that is extremely useful for confirming that the recalibration proceded as expected.  We highly recommend trying to rerun the script manually if possible.");
        }
    }
}

