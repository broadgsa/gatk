package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import Jama.Matrix;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.regex.Pattern;

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
    private final double[] empiricalMu; // BUGBUG: move these to VariantData class
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
        empiricalMu = new double[numAnnotations];
        empiricalSigma = new Matrix(numAnnotations, numAnnotations);
        isModelReadyForEvaluation = false;
    }

    public GaussianMixtureModel( final List<String> annotationLines, final List<String> gaussianLines ) { //BUGBUG: recreated here to match the integration tests

        gaussians = new ArrayList<MultivariateGaussian>( gaussianLines.size() );
        for( final String line : gaussianLines ) {
            final MultivariateGaussian gaussian = new MultivariateGaussian( annotationLines.size() );
            final String[] vals = line.split(",");
            gaussian.pMixtureLog10 = Math.log10( Double.parseDouble(vals[1]) );
            for( int iii = 0; iii < annotationLines.size(); iii++ ) {
                gaussian.mu[iii] = Double.parseDouble(vals[2+iii]); //BUGBUG: recreated here to match the integration tests
                for( int jjj = 0; jjj < annotationLines.size(); jjj++ ) {
                    gaussian.sigma.set(iii, jjj, Double.parseDouble(vals[2+annotationLines.size()+(iii*annotationLines.size())+jjj]) * 1.3); // BUGBUG: VRAC backOff, or get rid of this completely!?
                }
            }
            gaussians.add( gaussian );
        }

        this.shrinkage = 0.0; // not used when evaluating data, BUGBUG: move this to VariantData class
        this.dirichletParameter = 0.0; // not used when evaluating data
        empiricalMu = null; // not used when evaluating data
        empiricalSigma = null; // not used when evaluating data
        isModelReadyForEvaluation = false;
    }

    public void cacheEmpiricalStats( final List<VariantDatum> data ) {
        for( final VariantDatum datum : data ) {
            for( int jjj = 0; jjj < empiricalMu.length; jjj++ ) {
                empiricalMu[jjj] += datum.annotations[jjj] / ((double) data.size());
            }
        }
        empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, Matrix.identity(empiricalMu.length, empiricalMu.length)); // BUGBUG: why does the identity matrix work best here?
                               // is it because of a bug in the old implementation in which std>X variants were still counted in this calculation??
    }

    public void initializeRandomModel( final List<VariantDatum> data, final Random rand ) {

        // initialize random Gaussian means // BUGBUG: this is broken up this way to match the order of calls to rand.nextDouble() in the old code
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.initializeRandomMu( rand );
        }

        // initialize means using K-means algorithm
        initializeMeansUsingKMeans( data, 60, rand ); // BUGBUG: a VRAC argument?

        // initialize uniform mixture coefficients, random covariance matrices, and initial hyperparameters
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.pMixtureLog10 = Math.log10( 1.0 / ((double) gaussians.size()) );
            gaussian.initializeRandomSigma( rand );
            gaussian.hyperParameter_a = gaussian.mu.length;
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
                    if(dist < minDistance) {
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
                    if( datum.assignment == gaussian ) {
                        numAssigned++;
                        gaussian.incrementMu( datum );
                    }
                }
                if(numAssigned != 0) {
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
                likelihood += Math.pow( 10.0, pVarLog10 ); //BUGBUG: output recreated here to match the integration tests
            }
            final double[] pVarInGaussianReals = MathUtils.normalizeFromLog10( pVarInGaussianLog10 );
            int iii = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.assignPVarInGaussian( pVarInGaussianReals[iii++] ); //BUGBUG: to clean up
            }
        }

        logger.info("explained likelihood = " + String.format( "%.5f", likelihood / ((double) data.size()) ));
        return likelihood / ((double) data.size());
    }

    public void maximizationStep( final List<VariantDatum> data ) {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.maximizeGaussian( data, empiricalMu, empiricalSigma, shrinkage, dirichletParameter );
        }
    }

    private double getSumHyperParameterLambda() {
        double sum = 0.0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sum += gaussian.hyperParameter_lambda;
        }
        return sum;
    }

    public void output( final PrintStream clusterFile ) { //BUGBUG: output recreated here to match the integration tests
        normalizePMixtureLog10(); //BUGBUG: output recreated here to match the integration tests
        for( final MultivariateGaussian gaussian : gaussians ) {
            if( Math.pow(10.0, gaussian.pMixtureLog10) > 1E-4 ) {
                final double sigmaVals[][] = gaussian.sigma.getArray();
                clusterFile.print("@!CLUSTER");
                clusterFile.print(String.format(",%.8f", Math.pow(10.0, gaussian.pMixtureLog10)));
                for(int jjj = 0; jjj < gaussian.mu.length; jjj++ ) {
                    clusterFile.print(String.format(",%.8f", gaussian.mu[jjj]));
                }
                for(int jjj = 0; jjj < gaussian.mu.length; jjj++ ) {
                    for(int ppp = 0; ppp < gaussian.mu.length; ppp++ ) {
                        clusterFile.print(String.format(",%.8f", (sigmaVals[jjj][ppp] / gaussian.hyperParameter_a) )); // BUGBUG: this is a bug which should be fixed after passing integration tests
                    }
                }
                clusterFile.println();
            }
        }
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