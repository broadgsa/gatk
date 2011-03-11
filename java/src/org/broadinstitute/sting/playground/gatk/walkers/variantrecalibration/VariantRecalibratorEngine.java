package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.List;
import java.util.Random;

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

    private final static long RANDOM_SEED = 91801305;
    private final Random rand = new Random( RANDOM_SEED );
    private final static double MIN_PROB_CONVERGENCE = 1E-5;

    /////////////////////////////
    // Public Methods to interface with the Engine
    /////////////////////////////

    public VariantRecalibratorEngine( final VariantRecalibratorArgumentCollection VRAC ) {
        this.VRAC = VRAC;
        initialize( this.VRAC );
    }

    public GaussianMixtureModel generateModel( final List<VariantDatum> data ) {
        final GaussianMixtureModel model = new GaussianMixtureModel( 4, 3, 0.0001, 1000.0 ); //BUGBUG: VRAC.maxGaussians, VRAC.numAnnotations
        variationalBayesExpectationMaximization( model, data );
        return model;
    }

    public void evaluateData( final List<VariantDatum> data, final GaussianMixtureModel model ) {
        if( !model.isModelReadyForEvaluation ) {
            model.precomputeDenominatorForEvaluation();
        }
        for( final VariantDatum datum : data ) {
            datum.pVarGivenModel = evaluateDatum( datum, model );
        }
    }

    public void evaluateDataContrastively( final List<VariantDatum> data, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel ) {
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
        model.initializeRandomModel( data, rand );

        // The VBEM loop
        double previousLikelihood = -1E20;
        double currentLikelihood;
        int iteration = 1;
        while( iteration < 200 ) { //BUGBUG: VRAC.maxIterations
            currentLikelihood = model.expectationStep( data );
            model.maximizationStep( data );

            logger.info("Finished iteration " + iteration );
            iteration++;          
            if( Math.abs(currentLikelihood - previousLikelihood) < MIN_PROB_CONVERGENCE ) {
                logger.info("Convergence!");
                return; // Early return here because we have converged!
            }
            previousLikelihood = currentLikelihood;
        }
    }

    /////////////////////////////
    // Private Methods used for evaluating data given a GaussianMixtureModel
    /////////////////////////////

    private double evaluateDatum( final VariantDatum datum, final GaussianMixtureModel model ) {
        return model.evaluateDatum( datum );
    }
}
