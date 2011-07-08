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
import org.apache.commons.math.special.Gamma;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class MultivariateGaussian {
    public double pMixtureLog10;
    public double sumProb;
    final public double[] mu;
    final public Matrix sigma;
    public double hyperParameter_a;
    public double hyperParameter_b;
    public double hyperParameter_lambda;
    private double cachedDenomLog10;
    private Matrix cachedSigmaInverse;
    final private ExpandingArrayList<Double> pVarInGaussian;

    public MultivariateGaussian( final int numAnnotations ) {
        mu = new double[numAnnotations];
        sigma = new Matrix(numAnnotations, numAnnotations);
        pVarInGaussian = new ExpandingArrayList<Double>();
    }

    public void zeroOutMu() {
        Arrays.fill( mu, 0.0 );
    }

    public void zeroOutSigma() {
        final double[][] zeroSigma = new double[mu.length][mu.length];
        for( final double[] row : zeroSigma ) {
            Arrays.fill(row, 0);
        }
        final Matrix tmp = new Matrix(zeroSigma);
        sigma.setMatrix(0, mu.length - 1, 0, mu.length - 1, tmp);
    }

    public void initializeRandomMu( final Random rand ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] = -4.0 + 8.0 * rand.nextDouble();
        }
    }

    public void initializeRandomSigma( final Random rand ) {
        final double[][] randSigma = new double[mu.length][mu.length];
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = iii; jjj < mu.length; jjj++ ) {
                randSigma[jjj][iii] = 0.55 + 1.25 * rand.nextDouble();
                if( rand.nextBoolean() ) {
                    randSigma[jjj][iii] *= -1.0;
                }
                if( iii != jjj ) { randSigma[iii][jjj] = 0.0; } // Sigma is a symmetric, positive-definite matrix created by taking a lower diagonal matrix and multiplying it by its transpose
            }
        }
        Matrix tmp = new Matrix( randSigma );
        tmp = tmp.times(tmp.transpose());
        sigma.setMatrix(0, mu.length - 1, 0, mu.length - 1, tmp);
    }

    public double calculateDistanceFromMeanSquared( final VariantDatum datum ) {
        return MathUtils.distanceSquared( datum.annotations, mu );
    }

    public void incrementMu( final VariantDatum datum ) {
        incrementMu( datum, 1.0 );
    }
    
    public void incrementMu( final VariantDatum datum, final double prob ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] += prob * datum.annotations[jjj];
        }
    }

    public void divideEqualsMu( final double x ) {
        for( int jjj = 0; jjj < mu.length; jjj++ ) {
            mu[jjj] /= x;
        }
    }

    private void precomputeInverse() {
        try {
            cachedSigmaInverse = sigma.inverse();
        } catch( Exception e ) {
            throw new UserException("Error during clustering. Most likely there are too few variants used during Gaussian mixture modeling. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example).");
        }
    }


    public void precomputeDenominatorForEvaluation() {
        precomputeInverse();
        cachedDenomLog10 = Math.log10(Math.pow(2.0 * Math.PI, -1.0 * ((double) mu.length) / 2.0)) + Math.log10(Math.pow(sigma.det(), -0.5)) ;
    }

    public void precomputeDenominatorForVariationalBayes( final double sumHyperParameterLambda ) {

        // Variational Bayes calculations from Bishop
        precomputeInverse();
        cachedSigmaInverse.timesEquals( hyperParameter_a );
        double sum = 0.0;
        for(int jjj = 1; jjj <= mu.length; jjj++) {
            sum += Gamma.digamma( (hyperParameter_a + 1.0 - jjj) / 2.0 );
        }
        sum -= Math.log( sigma.det() );
        sum += Math.log(2.0) * mu.length;
        final double lambda = 0.5 * sum;
        final double pi = Gamma.digamma( hyperParameter_lambda ) - Gamma.digamma( sumHyperParameterLambda );
        final double beta = (-1.0 * mu.length) / (2.0 * hyperParameter_b);
        cachedDenomLog10 = (pi / Math.log(10.0)) + (lambda / Math.log(10.0)) + (beta / Math.log(10.0));
    }

    public double evaluateDatumLog10( final VariantDatum datum ) {
        double sumKernel = 0.0;
        final double[] crossProdTmp = new double[mu.length];
        Arrays.fill(crossProdTmp, 0.0);
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = 0; jjj < mu.length; jjj++ ) {
                crossProdTmp[iii] += (datum.annotations[jjj] - mu[jjj]) * cachedSigmaInverse.get(jjj, iii);
            }
        }
        for( int iii = 0; iii < mu.length; iii++ ) {
            sumKernel += crossProdTmp[iii] * (datum.annotations[iii] - mu[iii]);
        }
        
        return (( -0.5 * sumKernel ) / Math.log(10.0)) + cachedDenomLog10; // This is the definition of a Gaussian PDF Log10
    }

    public void assignPVarInGaussian( final double pVar ) {
        pVarInGaussian.add( pVar );
    }

    public void resetPVarInGaussian() {
        pVarInGaussian.clear();
    }

    public void maximizeGaussian( final List<VariantDatum> data, final double[] empiricalMu, final Matrix empiricalSigma,
                                  final double SHRINKAGE, final double DIRICHLET_PARAMETER, final double DEGREES_OF_FREEDOM ) {
        sumProb = 1E-10;
        final Matrix wishart = new Matrix(mu.length, mu.length);
        zeroOutMu();
        zeroOutSigma();
        
        int datumIndex = 0;
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            sumProb += prob;
            incrementMu( datum, prob );
        }
        divideEqualsMu( sumProb );

        final double shrinkageFactor = (SHRINKAGE * sumProb) / (SHRINKAGE + sumProb);
        for( int iii = 0; iii < mu.length; iii++ ) {
            for( int jjj = 0; jjj < mu.length; jjj++ ) {
                wishart.set(iii, jjj, shrinkageFactor * (mu[iii] - empiricalMu[iii]) * (mu[jjj] - empiricalMu[jjj]));
            }
        }

        datumIndex = 0;
        final Matrix pVarSigma = new Matrix(mu.length, mu.length);
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            for( int iii = 0; iii < mu.length; iii++ ) {
                for( int jjj = 0; jjj < mu.length; jjj++ ) {
                    pVarSigma.set(iii, jjj, prob * (datum.annotations[iii]-mu[iii]) * (datum.annotations[jjj]-mu[jjj]));
                }
            }
            sigma.plusEquals( pVarSigma );
        }

        sigma.plusEquals( empiricalSigma );
        sigma.plusEquals( wishart );

        for( int iii = 0; iii < mu.length; iii++ ) {
            mu[iii] = (sumProb * mu[iii] + SHRINKAGE * empiricalMu[iii]) / (sumProb + SHRINKAGE);
        }

        hyperParameter_a = sumProb + DEGREES_OF_FREEDOM;
        hyperParameter_b = sumProb + SHRINKAGE;
        hyperParameter_lambda = sumProb + DIRICHLET_PARAMETER;

        resetPVarInGaussian(); // clean up some memory
    }

    public void evaluateFinalModelParameters( final List<VariantDatum> data ) {
        sumProb = 0.0;
        zeroOutMu();
        zeroOutSigma();

        int datumIndex = 0;
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            sumProb += prob;
            incrementMu( datum, prob );
        }
        divideEqualsMu( sumProb );

        datumIndex = 0;
        final Matrix pVarSigma = new Matrix(mu.length, mu.length);
        for( final VariantDatum datum : data ) {
            final double prob = pVarInGaussian.get(datumIndex++);
            for( int iii = 0; iii < mu.length; iii++ ) {
                for( int jjj = 0; jjj < mu.length; jjj++ ) {
                    pVarSigma.set(iii, jjj, prob * (datum.annotations[iii]-mu[iii]) * (datum.annotations[jjj]-mu[jjj]));
                }
            }
            sigma.plusEquals( pVarSigma );
        }
        sigma.timesEquals( 1.0 / sumProb );

        resetPVarInGaussian(); // clean up some memory
    }
}