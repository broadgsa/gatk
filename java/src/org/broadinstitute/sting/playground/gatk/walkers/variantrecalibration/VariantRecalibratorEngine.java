package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantRecalibratorEngine {

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    protected final static Logger logger = Logger.getLogger(VariantRecalibratorEngine.class);

    // the unified argument collection
    final private VariantRecalibratorArgumentCollection VRAC;

    private final static double MIN_PROB_CONVERGENCE_LOG10 = 1.0;

    /////////////////////////////
    // Public Methods to interface with the Engine
    /////////////////////////////

    public VariantRecalibratorEngine( final VariantRecalibratorArgumentCollection VRAC ) {
        this.VRAC = VRAC;
        initialize( this.VRAC );
    }

    public GaussianMixtureModel generateModel( final List<VariantDatum> data ) {
        final GaussianMixtureModel model = new GaussianMixtureModel( VRAC.MAX_GAUSSIANS, data.get(0).annotations.length, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER );
        variationalBayesExpectationMaximization( model, data );
        return model;
    }

    public void evaluateData( final List<VariantDatum> data, final GaussianMixtureModel model, final boolean evaluateContrastively ) {
        if( !model.isModelReadyForEvaluation ) {
            model.precomputeDenominatorForEvaluation();
        }
        
        logger.info("Evaluating full set of " + data.size() + " variants...");
        for( final VariantDatum datum : data ) {
            final double thisLod = evaluateDatum( datum, model );
            datum.lod = ( evaluateContrastively ? (datum.prior + datum.lod - thisLod) : thisLod );
        }
    }

    /////////////////////////////
    // Private Methods used for initialization
    /////////////////////////////

    private void initialize( final VariantRecalibratorArgumentCollection VRAC ) {
    }

    /////////////////////////////
    // Private Methods used for generating a GaussianMixtureModel
    /////////////////////////////

    private void variationalBayesExpectationMaximization( final GaussianMixtureModel model, final List<VariantDatum> data ) {

        model.cacheEmpiricalStats( data );
        model.initializeRandomModel( data, VRAC.NUM_KMEANS_ITERATIONS );

        // The VBEM loop
        double previousLikelihood = model.expectationStep( data );
        double currentLikelihood;
        int iteration = 0;
        logger.info("Finished iteration " + iteration );
        while( iteration < VRAC.MAX_ITERATIONS ) {
            iteration++;
            model.maximizationStep( data );
            currentLikelihood = model.expectationStep( data );

            logger.info("Finished iteration " + iteration );
            if( Math.abs(currentLikelihood - previousLikelihood) < MIN_PROB_CONVERGENCE_LOG10 ) {
                logger.info("Convergence!");
                break;
            }
            previousLikelihood = currentLikelihood;
        }

        model.evaluateFinalModelParameters( data );
    }

    /////////////////////////////
    // Private Methods used for evaluating data given a GaussianMixtureModel
    /////////////////////////////

    private double evaluateDatum( final VariantDatum datum, final GaussianMixtureModel model ) {
        return model.evaluateDatum( datum );
    }
}
