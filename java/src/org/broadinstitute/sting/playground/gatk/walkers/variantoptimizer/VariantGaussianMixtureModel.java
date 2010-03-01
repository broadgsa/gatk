package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import java.util.Random;

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
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public class VariantGaussianMixtureModel extends VariantOptimizationModel {

    private final int numGaussians = 512;
    private final long RANDOM_SEED = 91801305;
    private final double MIN_PROB = 1E-30;
    private final double MIN_SUM_PROB = 1E-20;

    public VariantGaussianMixtureModel( VariantDataManager _dataManager, final double _targetTITV ) {
       super( _dataManager, _targetTITV );
    }
    
    public double[][] run() {

        final int numVariants = dataManager.numVariants;
        final int numAnnotations = dataManager.numAnnotations;
        final VariantDatum[] data = dataManager.data;
        final int numIterations = 16;

        final double[][] pTrueVariant = new double[numIterations][numVariants];
        final double[] pCluster = new double[numGaussians];
        final double[] probTi = new double[numGaussians];
        final double[] probTv = new double[numGaussians];
        final double[] clusterTruePositiveRate = new double[numGaussians];
        final double[][] pVarInCluster = new double[numGaussians][numVariants];
        final double[][] mu = new double[numGaussians][];
        final double[][] sigma = new double[numGaussians][];
        final Random rand = new Random( RANDOM_SEED );


        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            pCluster[kkk] = 1.0;
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            final double[] randSigma = new double[numAnnotations];
            if( dataManager.isNormalized ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = 0.9 + 0.2 * rand.nextDouble(); 
                }
            } else { // BUGBUG: if not normalized then the varianceVector hasn't been calculated --> null pointer
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = dataManager.varianceVector[jjj] + ((1.0 + rand.nextDouble()) * 0.01 * dataManager.varianceVector[jjj]);
                }
            }
            sigma[kkk] = randSigma;
        }

        for( int ttt = 0; ttt < numIterations; ttt++ ) {

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Expectation Step (calculate the probability that each data point is in each cluster)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for( int iii = 0; iii < numVariants; iii++ ) {
                double sumProb = 0.0;
                for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                    //pVarInCluster[kkk][iii] = rand.nextDouble();
                    double sum = 0.0;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        sum += ( (data[iii].annotations[jjj] - mu[kkk][jjj]) * (data[iii].annotations[jjj] - mu[kkk][jjj]) )
                                / sigma[kkk][jjj];
                    }
                    pVarInCluster[kkk][iii] = pCluster[kkk] * Math.exp( -0.5 * sum );

                    if( pVarInCluster[kkk][iii] < MIN_PROB) { // Very small numbers are a very big problem
                        pVarInCluster[kkk][iii] = MIN_PROB;
                    }

                    sumProb += pVarInCluster[kkk][iii];
                }

                if( sumProb > MIN_SUM_PROB ) { // Very small numbers are a very big problem
                    for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                        pVarInCluster[kkk][iii] /= sumProb;
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    mu[kkk][jjj] = 0.0;
                    sigma[kkk][jjj] = 0.0;
                }
            }
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                double sumProb = 0.0;
                for( int iii = 0; iii < numVariants; iii++ ) {
                    final double prob = pVarInCluster[kkk][iii];
                    sumProb += prob;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        mu[kkk][jjj] +=  prob * data[iii].annotations[jjj];
                    }
                }

                if( sumProb > MIN_SUM_PROB ) {
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        mu[kkk][jjj] /=  sumProb;
                    }
                }

                for( int iii = 0; iii < numVariants; iii++ ) {
                    final double prob = pVarInCluster[kkk][iii];
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        sigma[kkk][jjj] +=  prob * (data[iii].annotations[jjj]-mu[kkk][jjj]) * (data[iii].annotations[jjj]-mu[kkk][jjj]);
                    }
                }

                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    if( sigma[kkk][jjj] < MIN_PROB) { // Very small numbers are a very big problem
                        sigma[kkk][jjj] = MIN_PROB;
                    }
                }

                if( sumProb > MIN_SUM_PROB ) {
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        sigma[kkk][jjj] /=  sumProb;
                    }
                }

                pCluster[kkk] = sumProb / numVariants; // BUGBUG: Experiment with this, want to keep many clusters alive
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Evaluate the clusters using titv as an estimate of the true positive rate
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                probTi[kkk] = 0.0;
                probTv[kkk] = 0.0;
            }
            for( int iii = 0; iii < numVariants; iii++ ) {
                if( data[iii].isTransition ) { // transition
                    for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                        probTi[kkk] += pVarInCluster[kkk][iii];
                    }
                } else { // transversion
                    for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                        probTv[kkk] += pVarInCluster[kkk][iii];
                    }
                }
            }
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                clusterTruePositiveRate[kkk] = calcTruePositiveRateFromTITV( probTi[kkk] / probTv[kkk] );
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Evaluate each variant using the probability of being in each cluster and that cluster's true positive rate
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for( int iii = 0; iii < numVariants; iii++ ) {
                pTrueVariant[ttt][iii] = 0.0;
                for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                    pTrueVariant[ttt][iii] += pVarInCluster[kkk][iii] * clusterTruePositiveRate[kkk];
                }
            }



            System.out.println("Finished iteration " + (ttt+1) );
        }


        return pTrueVariant;
    }
}
