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

package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

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
import org.broadinstitute.sting.commandline.Argument;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing calibrated variant cluster parameters which can be applied to other datasets
 *
 * @author rpoplin
 * @since Feb 11, 2010
 *
 * @help.summary Takes variant calls as .vcf files, learns a Gaussian mixture model over the variant annotations producing calibrated variant cluster parameters which can be applied to other datasets
 */

public class VariantOptimizer extends RodWalker<ExpandingArrayList<VariantDatum>, ExpandingArrayList<VariantDatum>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="target_titv", shortName="titv", doc="The target Ti/Tv ratio towards which to optimize. (~~2.1 for whole genome experiments)", required=true)
    private double TARGET_TITV = 2.1;
    @Argument(fullName="ignore_all_input_filters", shortName="ignoreAllFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_ALL_INPUT_FILTERS = false;
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified the optimizer will use variants even if the specified filter name is marked in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="exclude_annotation", shortName="exclude", doc="The names of the annotations which should be excluded from the calculations", required=false)
    private String[] EXCLUDED_ANNOTATIONS = null;
    @Argument(fullName="force_annotation", shortName="force", doc="The names of the annotations which should be forced into the calculations even if they aren't present in every variant", required=false)
    private String[] FORCED_ANNOTATIONS = null;
    @Argument(fullName="clusterFile", shortName="clusterFile", doc="The output cluster file", required=true)
    private String CLUSTER_FILENAME = "optimizer.cluster";
    @Argument(fullName="numGaussians", shortName="nG", doc="The number of Gaussians to be used in the Gaussian Mixture model", required=false)
    private int NUM_GAUSSIANS = 7;
    @Argument(fullName="numIterations", shortName="nI", doc="The number of iterations to be performed in the Gaussian Mixture model", required=false)
    private int NUM_ITERATIONS = 10;
    @Argument(fullName="minVarInCluster", shortName="minVar", doc="The minimum number of variants in a cluster to be considered a valid cluster. It can be used to prevent overfitting.", required=false)
    private int MIN_VAR_IN_CLUSTER = 1000;
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";
    

    //@Argument(fullName="knn", shortName="knn", doc="The number of nearest neighbors to be used in the k-Nearest Neighbors model", required=false)
    //private int NUM_KNN = 2000;
    //@Argument(fullName = "optimization_model", shortName = "om", doc = "Optimization calculation model to employ -- GAUSSIAN_MIXTURE_MODEL is currently the default, while K_NEAREST_NEIGHBORS is also available for small callsets.", required = false)
    private VariantOptimizationModel.Model OPTIMIZATION_MODEL = VariantOptimizationModel.Model.GAUSSIAN_MIXTURE_MODEL;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final ExpandingArrayList<String> annotationKeys = new ExpandingArrayList<String>();
    private boolean firstVariant = true;
    private int numAnnotations = 0;
    private static final double INFINITE_ANNOTATION_VALUE = 10000.0;
    private Set<String> ignoreInputFilterSet = null;
    private boolean usingDBSNP = false;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet = new TreeSet<String>(Arrays.asList(IGNORE_INPUT_FILTERS));
        }
        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();
        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                usingDBSNP = true;
            }
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

        double annotationValues[] = new double[numAnnotations];

        for( final VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {

            if( vc != null && !vc.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) && vc.isSNP() ) {
                if( !vc.isFiltered() || IGNORE_ALL_INPUT_FILTERS || (ignoreInputFilterSet != null && ignoreInputFilterSet.containsAll(vc.getFilters())) ) {
                    if( firstVariant ) { // This is the first variant encountered so set up the list of annotations
                        annotationKeys.addAll( vc.getAttributes().keySet() );
                        annotationKeys.add("QUAL");
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
                        numAnnotations = annotationKeys.size();
                        annotationValues = new double[numAnnotations];
                        firstVariant = false;
                    }

                    int iii = 0;
                    for( final String key : annotationKeys ) {

                        double value = 0.0;
                        if( key.equals("AB") && !vc.getAttributes().containsKey(key) ) {
                            value = (0.5 - 0.005) + (0.01 * Math.random()); // HomVar calls don't have an allele balance
                        } else if( key.equals("QUAL") ) {
                            value = vc.getPhredScaledQual();
                        } else {
                            try {
                                value = Double.parseDouble( (String)vc.getAttribute( key, "0.0" ) );
                                if( Double.isInfinite(value) ) {
                                    value = ( value > 0 ? 1.0 : -1.0 ) * INFINITE_ANNOTATION_VALUE;
                                }
                            } catch( NumberFormatException e ) {
                                // do nothing, default value is 0.0
                            }
                        }
                        annotationValues[iii++] = value;
                    }

                    final VariantDatum variantDatum = new VariantDatum();
                    variantDatum.annotations = annotationValues;
                    variantDatum.isTransition = vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                    boolean isKnown = !vc.getAttribute("ID").equals(".");
                    if(usingDBSNP) {
                        isKnown = false;
                        for( VariantContext dbsnpVC : tracker.getVariantContexts(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, null, context.getLocation(), false, false) ) {
                            if(dbsnpVC != null && dbsnpVC.isSNP()) {
                                isKnown=true;
                            }
                        }
                    }
                    variantDatum.isKnown = isKnown;

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
                theModel = new VariantGaussianMixtureModel( dataManager, TARGET_TITV, NUM_GAUSSIANS, NUM_ITERATIONS, MIN_VAR_IN_CLUSTER );
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
            System.out.println( rScriptCommandLine );

            // Execute the RScript command to plot the table of truth values
            try {
                final Process p = Runtime.getRuntime().exec( rScriptCommandLine );
                p.waitFor();
            } catch (InterruptedException e) {
                e.printStackTrace();
                System.exit(-1);
            } catch ( IOException e ) {
                throw new StingException( "Unable to execute RScript command: " + rScriptCommandLine );
            }
        }
    }

}
