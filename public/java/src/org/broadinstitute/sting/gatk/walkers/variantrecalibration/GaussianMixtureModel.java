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

import Jama.Matrix;
import cern.jet.random.Normal;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class GaussianMixtureModel {

    protected final static Logger logger = Logger.getLogger(GaussianMixtureModel.class);

    private final ArrayList<MultivariateGaussian> gaussians;
    private final double shrinkage;
    private final double dirichletParameter;
    private final double priorCounts;
    private final double[] empiricalMu;
    private final Matrix empiricalSigma;
    public boolean isModelReadyForEvaluation;

    public GaussianMixtureModel( final int numGaussians, final int numAnnotations,
                                 final double shrinkage, final double dirichletParameter, final double priorCounts ) {

        gaussians = new ArrayList<MultivariateGaussian>( numGaussians );
        for( int iii = 0; iii < numGaussians; iii++ ) {
            final MultivariateGaussian gaussian = new MultivariateGaussian( numAnnotations );
            gaussians.add( gaussian );
        }
        this.shrinkage = shrinkage;
        this.dirichletParameter = dirichletParameter;
        this.priorCounts = priorCounts;
        empiricalMu = new double[numAnnotations];
        empiricalSigma = new Matrix(numAnnotations, numAnnotations);
        isModelReadyForEvaluation = false;
        Arrays.fill(empiricalMu, 0.0);
        empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, Matrix.identity(empiricalMu.length, empiricalMu.length).times(200.0).inverse());
    }

    public void initializeRandomModel( final List<VariantDatum> data, final int numKMeansIterations ) {

        // initialize random Gaussian means // BUGBUG: this is broken up this way to match the order of calls to rand.nextDouble() in the old code
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.initializeRandomMu( GenomeAnalysisEngine.getRandomGenerator() );
        }

        // initialize means using K-means algorithm
        logger.info( "Initializing model with " + numKMeansIterations + " k-means iterations..." );
        initializeMeansUsingKMeans( data, numKMeansIterations );

        // initialize uniform mixture coefficients, random covariance matrices, and initial hyperparameters
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.pMixtureLog10 = Math.log10( 1.0 / ((double) gaussians.size()) );
            gaussian.sumProb = 1.0 / ((double) gaussians.size());
            gaussian.initializeRandomSigma( GenomeAnalysisEngine.getRandomGenerator() );
            gaussian.hyperParameter_a = priorCounts;
            gaussian.hyperParameter_b = shrinkage;
            gaussian.hyperParameter_lambda = dirichletParameter;
        }
    }

    private void initializeMeansUsingKMeans( final List<VariantDatum> data, final int numIterations ) {

        int ttt = 0;
        while( ttt++ < numIterations ) {
            // E step: assign each variant to the nearest cluster
            for( final VariantDatum datum : data ) {
                double minDistance = Double.MAX_VALUE;
                MultivariateGaussian minGaussian = null;
                datum.assignment = minGaussian;
                for( final MultivariateGaussian gaussian : gaussians ) {
                    final double dist = gaussian.calculateDistanceFromMeanSquared( datum );
                    if( dist < minDistance ) {
                        minDistance = dist;
                        minGaussian = gaussian;
                    }
                }
                datum.assignment = minGaussian;
            }

            // M step: update gaussian means based on assigned variants
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.zeroOutMu();
                int numAssigned = 0;

                for( final VariantDatum datum : data ) {
                    if( datum.assignment.equals(gaussian) ) {
                        numAssigned++;
                        gaussian.incrementMu( datum );
                    }
                }
                if( numAssigned != 0 ) {
                    gaussian.divideEqualsMu( ((double) numAssigned) );
                } else {
                    gaussian.initializeRandomMu( GenomeAnalysisEngine.getRandomGenerator() );
                }
            }
        }
    }

    public void expectationStep( final List<VariantDatum> data ) {

        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForVariationalBayes( getSumHyperParameterLambda() );
        }

        for( final VariantDatum datum : data ) {
            final ArrayList<Double> pVarInGaussianLog10 = new ArrayList<Double>( gaussians.size() );
            for( final MultivariateGaussian gaussian : gaussians ) {
                final double pVarLog10 = gaussian.evaluateDatumLog10( datum );
                pVarInGaussianLog10.add( pVarLog10 );
            }
            final double[] pVarInGaussianNormalized = MathUtils.normalizeFromLog10( pVarInGaussianLog10 );
            int iii = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.assignPVarInGaussian( pVarInGaussianNormalized[iii++] ); //BUGBUG: to clean up
            }
        }
    }

    public void maximizationStep( final List<VariantDatum> data ) {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.maximizeGaussian( data, empiricalMu, empiricalSigma, shrinkage, dirichletParameter, priorCounts);
        }
    }

    private double getSumHyperParameterLambda() {
        double sum = 0.0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sum += gaussian.hyperParameter_lambda;
        }
        return sum;
    }

    public void evaluateFinalModelParameters( final List<VariantDatum> data ) {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.evaluateFinalModelParameters(data);
        }
        normalizePMixtureLog10();
    }

    public double normalizePMixtureLog10() {
        double sumDiff = 0.0;
        double sumPK = 0.0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sumPK += gaussian.sumProb;
        }

        int gaussianIndex = 0;
        double[] pGaussianLog10 = new double[gaussians.size()];
        for( final MultivariateGaussian gaussian : gaussians ) {
            pGaussianLog10[gaussianIndex++] = Math.log10( gaussian.sumProb / sumPK );
        }
        pGaussianLog10 = MathUtils.normalizeFromLog10( pGaussianLog10, true );

        gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sumDiff += Math.abs( pGaussianLog10[gaussianIndex] - gaussian.pMixtureLog10 );
            gaussian.pMixtureLog10 = pGaussianLog10[gaussianIndex++];
        }
        return sumDiff;
    }

    public void precomputeDenominatorForEvaluation() {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForEvaluation();
        }

        isModelReadyForEvaluation = true;
    }

    public double evaluateDatum( final VariantDatum datum ) {
        for( final boolean isNull : datum.isNull ) {
            if( isNull ) { return evaluateDatumMarginalized( datum ); }
        }
        // Fill an array with the log10 probability coming from each Gaussian and then use MathUtils to sum them up correctly
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + gaussian.evaluateDatumLog10( datum );
        }
        return MathUtils.log10sumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
    }

    // Used only to decide which covariate dimension is most divergent in order to report in the culprit info field annotation
    public Double evaluateDatumInOneDimension( final VariantDatum datum, final int iii ) {
        if(datum.isNull[iii]) { return null; }

        final Normal normal = new Normal(0.0, 1.0, null);
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            normal.setState( gaussian.mu[iii], gaussian.sigma.get(iii, iii) );
            pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + Math.log10( normal.pdf( datum.annotations[iii] ) );
        }
        return MathUtils.log10sumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
    }

    public double evaluateDatumMarginalized( final VariantDatum datum ) {
        int numRandomDraws = 0;
        double sumPVarInGaussian = 0.0;
        final int numIterPerMissingAnnotation = 10; // Trade off here between speed of computation and accuracy of the marginalization
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        // for each dimension
        for( int iii = 0; iii < datum.annotations.length; iii++ ) {
            // if it is missing marginalize over the missing dimension by drawing X random values for the missing annotation and averaging the lod
            if( datum.isNull[iii] ) {
                for( int ttt = 0; ttt < numIterPerMissingAnnotation; ttt++ ) {
                    datum.annotations[iii] = GenomeAnalysisEngine.getRandomGenerator().nextGaussian(); // draw a random sample from the standard normal distribution

                    // evaluate this random data point
                    int gaussianIndex = 0;
                    for( final MultivariateGaussian gaussian : gaussians ) {
                        pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + gaussian.evaluateDatumLog10( datum );
                    }

                    // add this sample's probability to the pile in order to take an average in the end
                    sumPVarInGaussian += Math.pow(10.0, MathUtils.log10sumLog10(pVarInGaussianLog10)); // p = 10 ^ Sum(pi_k * p(v|n,k))
                    numRandomDraws++;
                }
            }
        }
        return Math.log10( sumPVarInGaussian / ((double) numRandomDraws) );
    }
}