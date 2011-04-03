package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import Jama.Matrix;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

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
    private final double degreesOfFreedom;
    private final double[] empiricalMu;
    private final Matrix empiricalSigma;
    public boolean isModelReadyForEvaluation;

    public GaussianMixtureModel( final int numGaussians, final int numAnnotations,
                                 final double shrinkage, final double dirichletParameter ) {

        gaussians = new ArrayList<MultivariateGaussian>( numGaussians );
        for( int iii = 0; iii < numGaussians; iii++ ) {
            final MultivariateGaussian gaussian = new MultivariateGaussian( numAnnotations );
            gaussians.add( gaussian );
        }
        this.shrinkage = shrinkage;
        this.dirichletParameter = dirichletParameter;
        degreesOfFreedom = numAnnotations + 2;
        empiricalMu = new double[numAnnotations];
        empiricalSigma = new Matrix(numAnnotations, numAnnotations);
        isModelReadyForEvaluation = false;
    }

    public void cacheEmpiricalStats( final List<VariantDatum> data ) {
        final double[][] tmpSigmaVals = new double[empiricalMu.length][empiricalMu.length];
        for( int iii = 0; iii < empiricalMu.length; iii++ ) {
            empiricalMu[iii] = 0.0;
            for( int jjj = iii; jjj < empiricalMu.length; jjj++ ) {
                tmpSigmaVals[iii][jjj] = 0.0;
            }
        }

        for( final VariantDatum datum : data ) {
            for( int iii = 0; iii < empiricalMu.length; iii++ ) {
                empiricalMu[iii] += datum.annotations[iii] / ((double) data.size());
            }
        }

        //for( final VariantDatum datum : data ) {
        //    for( int iii = 0; iii < empiricalMu.length; iii++ ) {
        //        for( int jjj = 0; jjj < empiricalMu.length; jjj++ ) {
        //            tmpSigmaVals[iii][jjj] += (datum.annotations[iii]-empiricalMu[iii]) * (datum.annotations[jjj]-empiricalMu[jjj]);
        //        }
        //    }
        //}

        //empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, new Matrix(tmpSigmaVals));
        //empiricalSigma.timesEquals( 1.0 / ((double) data.size()) );
        //empiricalSigma.timesEquals( 1.0 / (Math.pow(gaussians.size(), 2.0 / ((double) empiricalMu.length))) );
        empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, Matrix.identity(empiricalMu.length, empiricalMu.length));
    }

    public void initializeRandomModel( final List<VariantDatum> data, final Random rand, final int numKMeansIterations ) {

        // initialize random Gaussian means // BUGBUG: this is broken up this way to match the order of calls to rand.nextDouble() in the old code
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.initializeRandomMu( rand );
        }

        // initialize means using K-means algorithm
        logger.info( "Initializing model with " + numKMeansIterations + " k-means iterations..." );
        initializeMeansUsingKMeans( data, numKMeansIterations, rand );

        // initialize uniform mixture coefficients, random covariance matrices, and initial hyperparameters
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.pMixtureLog10 = Math.log10( 1.0 / ((double) gaussians.size()) );
            gaussian.initializeRandomSigma( rand );
            gaussian.hyperParameter_a = degreesOfFreedom;
            gaussian.hyperParameter_b = shrinkage;
            gaussian.hyperParameter_lambda = dirichletParameter;
        }
    }

    private void initializeMeansUsingKMeans( final List<VariantDatum> data, final int numIterations, final Random rand ) {

        int ttt = 0;
        while( ttt++ < numIterations ) {
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
                    gaussian.initializeRandomMu( rand );
                }
            }
        }
    }

    public double expectationStep( final List<VariantDatum> data ) {

        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForVariationalBayes( getSumHyperParameterLambda() );
        }

        double likelihood = 0.0;
        for( final VariantDatum datum : data ) {
            final ArrayList<Double> pVarInGaussianLog10 = new ArrayList<Double>( gaussians.size() );
            for( final MultivariateGaussian gaussian : gaussians ) {
                final double pVarLog10 = gaussian.evaluateDatumLog10( datum );
                pVarInGaussianLog10.add( pVarLog10 );
                likelihood += pVarLog10;
            }
            final double[] pVarInGaussianNormalized = MathUtils.normalizeFromLog10( pVarInGaussianLog10 );
            int iii = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.assignPVarInGaussian( pVarInGaussianNormalized[iii++] ); //BUGBUG: to clean up
            }
        }

        final double scaledTotalLikelihoodLog10 = likelihood / ((double) data.size());
        logger.info( "sum Log10 likelihood = " + String.format("%.5f", scaledTotalLikelihoodLog10) );
        return scaledTotalLikelihoodLog10;
    }

    public void maximizationStep( final List<VariantDatum> data ) {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.maximizeGaussian( data, empiricalMu, empiricalSigma, shrinkage, dirichletParameter, degreesOfFreedom );
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
            gaussian.evaluateFinalModelParameters( data );
        }
        normalizePMixtureLog10();
    }

    private void normalizePMixtureLog10() {
        double sumPK = 0.0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sumPK += gaussian.pMixtureLog10;
        }

        int gaussianIndex = 0;
        double[] pGaussianLog10 = new double[gaussians.size()];
        for( final MultivariateGaussian gaussian : gaussians ) {
            pGaussianLog10[gaussianIndex++] = Math.log10( gaussian.pMixtureLog10 / sumPK ); //BUGBUG: to clean up
        }
        pGaussianLog10 = MathUtils.normalizeFromLog10( pGaussianLog10, true );

        gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.pMixtureLog10 = pGaussianLog10[gaussianIndex++];
        }
    }

    public void precomputeDenominatorForEvaluation() {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForEvaluation();
        }

        isModelReadyForEvaluation = true;
    }

    public double evaluateDatum( final VariantDatum datum ) {
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + gaussian.evaluateDatumLog10( datum );
        }
        return MathUtils.log10sumLog10(pVarInGaussianLog10);
    }
}