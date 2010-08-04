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
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.commandline.Argument;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing calibrated variant cluster parameters which can be applied to other datasets
 *
 * @author rpoplin
 * @since Feb 11, 2010
 *
 * @help.summary Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing calibrated variant cluster parameters which can be applied to other datasets
 */

public class GenerateVariantClustersWalker extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="ignore_all_input_filters", shortName="ignoreAllFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_ALL_INPUT_FILTERS = false;
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the optimizer will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="use_annotation", shortName="an", doc="The names of the annotations which should used for calculations", required=true)
    private String[] USE_ANNOTATIONS = null;
    @Argument(fullName="clusterFile", shortName="clusterFile", doc="The output cluster file", required=true)
    private String CLUSTER_FILENAME = "optimizer.cluster";
    @Argument(fullName="maxGaussians", shortName="mG", doc="The maximum number of Gaussians to try during Bayesian clustering", required=false)
    private int MAX_GAUSSIANS = 30;
    @Argument(fullName="maxIterations", shortName="mI", doc="The maximum number of iterations to be performed when clustering. Clustering will normally end when convergence is detected.", required=false)
    private int MAX_ITERATIONS = 200;
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is maybe /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";
    @Argument(fullName="weightNovels", shortName="weightNovels", doc="The weight for novel variants during clustering", required=false)
    private double WEIGHT_NOVELS = 0.0;
    @Argument(fullName="weightKnowns", shortName="weightKnowns", doc="The weight for MQ2+ known variants during clustering", required=false)
    private double WEIGHT_KNOWNS = 0.0;
    @Argument(fullName="weightHapMap", shortName="weightHapMap", doc="The weight for known HapMap variants during clustering", required=false)
    private double WEIGHT_HAPMAP = 1.0;
    @Argument(fullName="weight1000Genomes", shortName="weight1000Genomes", doc="The weight for known 1000 Genomes Project variants during clustering", required=false)
    private double WEIGHT_1000GENOMES = 1.0;
    @Argument(fullName="weightMQ1", shortName="weightMQ1", doc="The weight for MQ1 dbSNP variants during clustering", required=false)
    private double WEIGHT_MQ1 = 0.0;
    @Argument(fullName="forceIndependent", shortName="forceIndependent", doc="Force off-diagonal entries in the covariance matrix to be zero.", required=false)
    private boolean FORCE_INDEPENDENT = false;
    @Argument(fullName="stdThreshold", shortName="std", doc="If a variant has annotations more than -std standard deviations away from mean then don't use it for clustering.", required=false)
    private double STD_THRESHOLD = 6.0;
    @Argument(fullName="qualThreshold", shortName="qual", doc="If a known variant has raw QUAL value less than -qual then don't use it for clustering.", required=false)
    private double QUAL_THRESHOLD = 300.0;
    @Argument(fullName="shrinkage", shortName="shrinkage", doc="The shrinkage parameter in variational Bayes algorithm.", required=false)
    private double SHRINKAGE = 0.0001;
    @Argument(fullName="dirichlet", shortName="dirichlet", doc="The dirichlet parameter in variational Bayes algoirthm.", required=false)
    private double DIRICHLET_PARAMETER = 1000.0;


    //@Argument(fullName="knn", shortName="knn", doc="The number of nearest neighbors to be used in the k-Nearest Neighbors model", required=false)
    //private int NUM_KNN = 2000;
    //@Argument(fullName = "optimization_model", shortName = "om", doc = "Optimization calculation model to employ -- GAUSSIAN_MIXTURE_MODEL is currently the default, while K_NEAREST_NEIGHBORS is also available for small callsets.", required = false)
    private VariantOptimizationModel.Model OPTIMIZATION_MODEL = VariantOptimizationModel.Model.GAUSSIAN_MIXTURE_MODEL;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private ExpandingArrayList<String> annotationKeys;
    private Set<String> ignoreInputFilterSet = null;
    private int maxAC = 0;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }
        
        annotationKeys = new ExpandingArrayList<String>(Arrays.asList(USE_ANNOTATIONS));

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet = new TreeSet<String>(Arrays.asList(IGNORE_INPUT_FILTERS));
        }

        boolean foundDBSNP = false;
        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();
        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                foundDBSNP = true;
            }
        }

        if(!foundDBSNP) {
            throw new StingException("dbSNP track is required. This calculation is critically dependent on being able to distinguish known and novel sites.");
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

        final double annotationValues[] = new double[annotationKeys.size()];

        // todo -- do we really need to support multiple tracks -- logic is cleaner without this case -- what's the use case?
        for( final VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {
            if( vc != null && !vc.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) && vc.isSNP() ) {
                if( !vc.isFiltered() || IGNORE_ALL_INPUT_FILTERS || (ignoreInputFilterSet != null && ignoreInputFilterSet.containsAll(vc.getFilters())) ) {
                    int iii = 0;
                    for( final String key : annotationKeys ) {
                        annotationValues[iii++] = VariantGaussianMixtureModel.decodeAnnotation( key, vc, true );
                    }

                    final VariantDatum variantDatum = new VariantDatum();
                    variantDatum.annotations = annotationValues;
                    variantDatum.isTransition = vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                    variantDatum.alleleCount = vc.getChromosomeCount(vc.getAlternateAllele(0)); // BUGBUG: assumes file has genotypes
                    if( variantDatum.alleleCount > maxAC ) {
                        maxAC = variantDatum.alleleCount;
                    }

                    variantDatum.isKnown = false;
                    variantDatum.weight = WEIGHT_NOVELS;
                    variantDatum.qual = vc.getPhredScaledQual();

                    final DbSNPFeature dbsnp = DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
                    if( dbsnp != null ) {
                        variantDatum.isKnown = true;
                        variantDatum.weight = WEIGHT_KNOWNS;

                        if( DbSNPHelper.isHapmap( dbsnp ) ) { variantDatum.weight = WEIGHT_HAPMAP; }
                        else if( DbSNPHelper.is1000genomes( dbsnp ) ) { variantDatum.weight = WEIGHT_1000GENOMES; }
                        else if( DbSNPHelper.isMQ1( dbsnp ) ) { variantDatum.weight = WEIGHT_MQ1; }
                    }
                    
                    mapList.add( variantDatum );
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

        final VariantDataManager dataManager = new VariantDataManager( reduceSum, annotationKeys );
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        logger.info( "There are " + dataManager.numVariants + " variants and " + dataManager.numAnnotations + " annotations." );
        logger.info( "The annotations are: " + annotationKeys );

        dataManager.normalizeData(); // Each data point is now [ (x - mean) / standard deviation ]

        // Create either the Gaussian Mixture Model or the Nearest Neighbors model and run it
        VariantGaussianMixtureModel theModel;
        switch (OPTIMIZATION_MODEL) {
            case GAUSSIAN_MIXTURE_MODEL:
                theModel = new VariantGaussianMixtureModel( dataManager, MAX_GAUSSIANS, MAX_ITERATIONS, maxAC, FORCE_INDEPENDENT,
                                                            STD_THRESHOLD, QUAL_THRESHOLD, SHRINKAGE, DIRICHLET_PARAMETER );
                break;
            //case K_NEAREST_NEIGHBORS:
            //    theModel = new VariantNearestNeighborsModel( dataManager, TARGET_TITV, NUM_KNN );
            //    break;
            default:
                throw new StingException( "Variant Optimization Model is unrecognized. Implemented options are GAUSSIAN_MIXTURE_MODEL and K_NEAREST_NEIGHBORS" );
        }
        
        theModel.run( CLUSTER_FILENAME );
        theModel.outputClusterReports( CLUSTER_FILENAME );

        for( final String annotation : annotationKeys ) {
            // Execute Rscript command to plot the optimization curve
            // Print out the command line to make it clear to the user what is being executed and how one might modify it
            final String rScriptCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_ClusterReport.R" + " " + CLUSTER_FILENAME + "." + annotation + ".dat " + annotation;
            logger.info( rScriptCommandLine );

            // Execute the RScript command to plot the table of truth values
            try {
                Runtime.getRuntime().exec( rScriptCommandLine );
            } catch ( IOException e ) {
                Utils.warnUser("Unable to execute the RScript command because of [" + e.getMessage() + "].  While not critical to the calculations themselves, the script outputs a report that is extremely useful for confirming that the clustering proceded as expected.  We highly recommend trying to rerun the script manually if possible.");
            }
        }
    }

}
