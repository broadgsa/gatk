/*
 * Copyright (c) 2011 The Broad Institute
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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;

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
    public final static double MIN_ACCEPTABLE_LOD_SCORE = -20000.0;

    // the unified argument collection
    final private VariantRecalibratorArgumentCollection VRAC;

    private final static double MIN_PROB_CONVERGENCE = 2E-2;

    /////////////////////////////
    // Public Methods to interface with the Engine
    /////////////////////////////

    public VariantRecalibratorEngine( final VariantRecalibratorArgumentCollection VRAC ) {
        this.VRAC = VRAC;
    }

    public GaussianMixtureModel generateModel( final List<VariantDatum> data ) {
        final GaussianMixtureModel model = new GaussianMixtureModel( VRAC.MAX_GAUSSIANS, data.get(0).annotations.length, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS );
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
            if( Double.isNaN(thisLod) ) {
                throw new UserException("NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example)");
            }

            datum.lod = ( evaluateContrastively ?
                            ( Double.isInfinite(datum.lod) ? // positive model said negative infinity
                                    ( MIN_ACCEPTABLE_LOD_SCORE + GenomeAnalysisEngine.getRandomGenerator().nextDouble() * MIN_ACCEPTABLE_LOD_SCORE ) // Negative infinity lod values are possible when covariates are extremely far away from their tight Gaussians
                                    : datum.prior + datum.lod - thisLod) // contrastive evaluation: (prior + positive model - negative model)
                            : thisLod ); // positive model only so set the lod and return
        }
    }

    public void calculateWorstPerformingAnnotation( final List<VariantDatum> data, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel ) {
        for( final VariantDatum datum : data ) {
            int worstAnnotation = -1;
            double minProb = Double.MAX_VALUE;
            for( int iii = 0; iii < datum.annotations.length; iii++ ) {
                final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
                final Double badProbLog10 = badModel.evaluateDatumInOneDimension(datum, iii);
                if( goodProbLog10 != null && badProbLog10 != null ) {
                    final double prob = goodProbLog10 - badProbLog10;
                    if(prob < minProb) { minProb = prob; worstAnnotation = iii; }
                }
            }
            datum.worstAnnotation = worstAnnotation;
        }
    }


    /////////////////////////////
    // Private Methods used for generating a GaussianMixtureModel
    /////////////////////////////

    private void variationalBayesExpectationMaximization( final GaussianMixtureModel model, final List<VariantDatum> data ) {

        model.initializeRandomModel( data, VRAC.NUM_KMEANS_ITERATIONS );

        // The VBEM loop
        model.normalizePMixtureLog10();
        model.expectationStep( data );
        double currentChangeInMixtureCoefficients;
        int iteration = 0;
        logger.info("Finished iteration " + iteration + ".");
        while( iteration < VRAC.MAX_ITERATIONS ) {
            iteration++;
            model.maximizationStep( data );
            currentChangeInMixtureCoefficients = model.normalizePMixtureLog10();
            model.expectationStep(data);
            if( iteration % 5 == 0 ) { // cut down on the number of output lines so that users can read the warning messages
                logger.info("Finished iteration " + iteration + ". \tCurrent change in mixture coefficients = " + String.format("%.5f", currentChangeInMixtureCoefficients));
            }
            if( iteration > 2 && currentChangeInMixtureCoefficients < MIN_PROB_CONVERGENCE ) {
                logger.info("Convergence after " + iteration + " iterations!");
                break;
            }
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
