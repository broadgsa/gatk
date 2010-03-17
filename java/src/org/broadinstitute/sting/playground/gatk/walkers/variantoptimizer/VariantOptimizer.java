package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.IOException;

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
    @Argument(fullName="target_titv", shortName="titv", doc="The target Ti/Tv ratio towards which to optimize. (~~2.2 for whole genome experiments)", required=true)
    private double TARGET_TITV = 2.1;
    @Argument(fullName="ignore_input_filters", shortName="ignoreFilters", doc="If specified the optimizer will use variants even if the FILTER column is marked in the VCF file", required=false)
    private boolean IGNORE_INPUT_FILTERS = false;
    @Argument(fullName="exclude_annotation", shortName="exclude", doc="The names of the annotations which should be excluded from the calculations", required=false)
    private String[] EXCLUDED_ANNOTATIONS = null;
    @Argument(fullName="force_annotation", shortName="force", doc="The names of the annotations which should be forced into the calculations even if they aren't present in every variant", required=false)
    private String[] FORCED_ANNOTATIONS = null;
    @Argument(fullName="clusterFile", shortName="clusterFile", doc="The output cluster file", required=true)
    private String CLUSTER_FILENAME = "optimizer.cluster";
    @Argument(fullName="numGaussians", shortName="nG", doc="The number of Gaussians to be used in the Gaussian Mixture model", required=false)
    private int NUM_GAUSSIANS = 32;
    @Argument(fullName="numIterations", shortName="nI", doc="The number of iterations to be performed in the Gaussian Mixture model", required=false)
    private int NUM_ITERATIONS = 10;
    @Argument(fullName="minVarInCluster", shortName="minVar", doc="The minimum number of variants in a cluster to be considered a valid cluster. Used to prevent overfitting. Default is 2000.", required=true)
    private int MIN_VAR_IN_CLUSTER = 2000;

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

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        //if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }   
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
                        // do nothing, default value is 0.0
                    }
                    annotationValues[iii++] = value;
                }

                // Variant quality ("QUAL") is not in the list of annotations, but is useful so add it here.
                annotationValues[iii] = vc.getPhredScaledQual();

                VariantDatum variantDatum = new VariantDatum();
                variantDatum.annotations = annotationValues;
                variantDatum.isTransition = vc.getSNPSubstitutionType().compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0;
                variantDatum.isKnown = !vc.getAttribute("ID").equals(".");
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

        final VariantDataManager dataManager = new VariantDataManager( reduceSum, annotationKeys );
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        logger.info( "There are " + dataManager.numVariants + " variants and " + dataManager.numAnnotations + " annotations." );
        logger.info( "The annotations are: " + annotationKeys + " and QUAL." );

        dataManager.normalizeData(); // Each data point is now [ (x - mean) / standard deviation ]
        
        // Create either the Gaussian Mixture Model or the Nearest Neighbors model and run it
        VariantOptimizationModel theModel;
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
    }

}
