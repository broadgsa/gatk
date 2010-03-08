package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import java.io.PrintStream;
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

public final class VariantGaussianMixtureModel extends VariantOptimizationModel implements VariantClusteringModel {

    private final int numGaussians;
    private final int numIterations;
    private final long RANDOM_SEED = 91801305;
    private final Random rand = new Random( RANDOM_SEED );
    private final double MIN_PROB = 1E-30;
    private final double MIN_SUM_PROB = 1E-20;

    private final double[][] mu;
    private final double[][] sigma;
    private final double[] pCluster;
    private final int[] numMaxClusterKnown;
    private final int[] numMaxClusterNovel;
    private final double[] clusterTITV;
    private final double[] clusterTruePositiveRate;

    public VariantGaussianMixtureModel( VariantDataManager _dataManager, final double _targetTITV, final int _numGaussians, final int _numIterations ) {
       super( _dataManager, _targetTITV );
        numGaussians = _numGaussians;
        numIterations = _numIterations;

        mu = new double[numGaussians][];
        sigma = new double[numGaussians][];
        pCluster = new double[numGaussians];
        numMaxClusterKnown = new int[numGaussians];
        numMaxClusterNovel = new int[numGaussians];
        clusterTITV = new double[numGaussians];
        clusterTruePositiveRate = new double[numGaussians];
    }
    
    public final void run( final String outputPrefix ) {

        // Create the subset of the data to cluster with
        int numNovel = 0;
        for( final VariantDatum datum : dataManager.data ) {
            if( !datum.isKnown ) {
                numNovel++;
            }
        }
        VariantDatum[] data;

        if( numNovel * 2 * 1.3 < dataManager.numVariants ) {
            data = new VariantDatum[numNovel*2];
            int iii = 0;
            for( final VariantDatum datum : dataManager.data ) {
                if( !datum.isKnown ) {
                    data[iii++] = datum;
                }
            }
            while( iii < numNovel*2 ) { // grab an equal number of known variants at random
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isKnown ) {
                    data[iii++] = datum;
                }
            }
        } else {
            data = dataManager.data;
        }

        System.out.println("Clustering with " + numNovel + " novel variants and " + (data.length - numNovel) + " known variants...");
        if( data.length == dataManager.numVariants ) { System.out.println(" (used all variants since 2*numNovel is so large compared to the full set) "); }
        createClusters( data ); // Using a subset of the data
        System.out.println("Printing out cluster parameters...");
        printClusters( outputPrefix );
        System.out.println("Applying clusters to all the variants...");
        applyClusters( dataManager.data, outputPrefix ); // Using all the data
    }

    public final void createClusters( final VariantDatum[] data ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[][] pVarInCluster = new double[numGaussians][numVariants];
        final double[] probTi = new double[numGaussians];
        final double[] probTv = new double[numGaussians];

        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        // Set up the initial random Gaussians
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            numMaxClusterKnown[kkk] = 0;
            numMaxClusterNovel[kkk] = 0;
            pCluster[kkk] = 1.0 / ((double) numGaussians);
            //final double[] randMu = new double[numAnnotations];
            //for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            //    randMu[jjj] = data[rand.nextInt(numVariants)].annotations[jjj];
            //}
            mu[kkk] = data[rand.nextInt(numVariants)].annotations; //randMu;
            final double[] randSigma = new double[numAnnotations];
            if( dataManager.isNormalized ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = 0.75 + 0.4 * rand.nextDouble();
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
            evaluateGaussians( data, pVarInCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            maximizeGaussians( data, pVarInCluster );

            System.out.println("Finished iteration " + (ttt+1) );
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Evaluate the clusters using titv as an estimate of the true positive rate
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        evaluateGaussians( data, pVarInCluster ); // One final evaluation because the Gaussians moved in the last maximization step

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

            double maxProb = pVarInCluster[0][iii];
            int maxCluster = 0;
            for( int kkk = 1; kkk < numGaussians; kkk++ ) {
                if( pVarInCluster[kkk][iii] > maxProb ) {
                    maxProb = pVarInCluster[kkk][iii];
                    maxCluster = kkk;
                }
            }
            if( data[iii].isKnown ) {
                numMaxClusterKnown[maxCluster]++;
            } else {
                numMaxClusterNovel[maxCluster]++;
            }
        }
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            clusterTITV[kkk] = probTi[kkk] / probTv[kkk];
            clusterTruePositiveRate[kkk] = calcTruePositiveRateFromTITV( clusterTITV[kkk] );
        }
    }

    private void printClusters( final String outputPrefix ) {
        try {
            final PrintStream outputFile = new PrintStream( outputPrefix + ".clusters" );
            int clusterNumber = 0;
            final int numAnnotations = mu[0].length;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( numMaxClusterKnown[kkk] + numMaxClusterNovel[kkk] >= 2000 ) {
                    outputFile.print(clusterNumber + ",");
                    outputFile.print(numMaxClusterKnown[kkk] + ",");
                    outputFile.print(numMaxClusterNovel[kkk] + ",");
                    outputFile.print(clusterTITV[kkk] + ",");
                    outputFile.print(clusterTruePositiveRate[kkk] + ",");
                    for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        outputFile.print(mu[kkk][jjj] + ",");
                    }
                    for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        outputFile.print(sigma[kkk][jjj] + ",");
                    }
                    outputFile.println(-1);
                    clusterNumber++;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public final void applyClusters( final VariantDatum[] data, final String outputPrefix ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[] pTrueVariant = new double[numVariants];
        final double[] pVarInCluster = new double[numGaussians];
        final int[] clusterAssignment = new int[numVariants];

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Evaluate each variant using the probability of being in each cluster and that cluster's true positive rate
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for( int iii = 0; iii < numVariants; iii++ ) {
            evaluateGaussiansForSingleVariant( data[iii], pVarInCluster );

            pTrueVariant[iii] = 0.0;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                pTrueVariant[iii] += pVarInCluster[kkk] * clusterTruePositiveRate[kkk];
            }

            double maxProb = -1.0;//pVarInCluster[0][iii];
            int maxCluster = -1;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( numMaxClusterKnown[kkk] + numMaxClusterNovel[kkk] >= 2000 && pVarInCluster[kkk] > maxProb ) {
                    maxProb = pVarInCluster[kkk];
                    maxCluster = kkk;
                }
            }
            clusterAssignment[iii] = maxCluster;
        }

        PrintStream outputFile = null;
        try {
            outputFile = new PrintStream( outputPrefix + ".data" );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        for(int iii = 0; iii < numVariants; iii++) {
            outputFile.print(pTrueVariant[iii] + ",");
            outputFile.print(clusterAssignment[iii] + ",");
            for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                outputFile.print(data[iii].originalAnnotations[jjj] + ",");
            }
            outputFile.println( (data[iii].isTransition ? 1 : 0)
                    + "," + (data[iii].isKnown? 1 : 0)
                    + "," + (data[iii].isFiltered ? 1 : 0));
        }

        try {
            outputFile = new PrintStream( outputPrefix + ".optimize" );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        for(int iii = 0; iii < numVariants; iii++) {
            outputFile.print(String.format("%.4f",pTrueVariant[iii]) + ",");
            outputFile.println( (data[iii].isTransition ? 1 : 0)
                    + "," + (data[iii].isKnown? 1 : 0)
                    + "," + (data[iii].isFiltered ? 1 : 0));
        }
    }


    private void evaluateGaussians( final VariantDatum[] data, final double[][] pVarInCluster ) {

        final int numAnnotations = data[0].annotations.length;

        for( int iii = 0; iii < data.length; iii++ ) {
            double sumProb = 0.0;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
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
    }


    private void evaluateGaussiansForSingleVariant( final VariantDatum datum, final double[] pVarInCluster ) {

        final int numAnnotations = datum.annotations.length;

        double sumProb = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            double sum = 0.0;
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                sum += ( (datum.annotations[jjj] - mu[kkk][jjj]) * (datum.annotations[jjj] - mu[kkk][jjj]) )
                        / sigma[kkk][jjj];
            }
            pVarInCluster[kkk] = pCluster[kkk] * Math.exp( -0.5 * sum );

            if( pVarInCluster[kkk] < MIN_PROB) { // Very small numbers are a very big problem
                pVarInCluster[kkk] = MIN_PROB;
            }

            sumProb += pVarInCluster[kkk];
        }

        if( sumProb > MIN_SUM_PROB ) { // Very small numbers are a very big problem
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                pVarInCluster[kkk] /= sumProb;
            }
        }
    }


    private void maximizeGaussians( final VariantDatum[] data, final double[][] pVarInCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

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

            pCluster[kkk] = sumProb / numVariants;
        }

        // Clean up extra big or extra small clusters
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            if( pCluster[kkk] > 0.45 ) { // This is a very large cluster compared to all the others
                System.out.println("Very large cluster!");
                final int numToReplace = 4;
                final double[] savedSigma = sigma[kkk];
                for( int rrr = 0; rrr < numToReplace; rrr++ ) {
                    // Find an example variant in the large cluster, drawn randomly
                    int randVarIndex = -1;
                    boolean foundVar = false;
                    while( !foundVar ) {
                        randVarIndex = rand.nextInt( numVariants );
                        final double probK = pVarInCluster[kkk][randVarIndex];
                        boolean inClusterK = true;
                        for( int ccc = 0; ccc < numGaussians; ccc++ ) {
                            if( pVarInCluster[ccc][randVarIndex] > probK ) {
                                inClusterK = false;
                                break;
                            }
                        }
                        if( inClusterK ) { foundVar = true; }
                    }

                    // Find a place to put the example variant
                    if( rrr == 0 ) { // Replace the big cluster that kicked this process off
                        mu[kkk] = data[randVarIndex].annotations;
                        pCluster[kkk] = 1.0 / ((double) numGaussians);
                    } else { // Replace the cluster with the minimum prob
                        double minProb = pCluster[0];
                        int minClusterIndex = 0;
                        for( int ccc = 1; ccc < numGaussians; ccc++ ) {
                            if( pCluster[ccc] < minProb ) {
                                minProb = pCluster[ccc];
                                minClusterIndex = ccc;
                            }
                        }
                        mu[minClusterIndex] = data[randVarIndex].annotations;
                        sigma[minClusterIndex] = savedSigma;
                        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                            sigma[minClusterIndex][jjj] += -0.06 + 0.12 * rand.nextDouble();
                            if( sigma[minClusterIndex][jjj] < MIN_SUM_PROB ) {
                                sigma[minClusterIndex][jjj] = MIN_SUM_PROB;
                            }
                        }
                        pCluster[minClusterIndex] = 1.0 / ((double) numGaussians);
                    }
                }
            }
        }

        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            if( pCluster[kkk] < 0.05 * (1.0 / ((double) numGaussians)) ) { // This is a very small cluster compared to all the others
                System.out.println("Very small cluster!");
                pCluster[kkk] = 1.0 / ((double) numGaussians);
                mu[kkk] = data[rand.nextInt(numVariants)].annotations;
                final double[] randSigma = new double[numAnnotations];
                if( dataManager.isNormalized ) {
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        randSigma[jjj] = 0.6 + 0.4 * rand.nextDouble(); // Explore a wider range
                    }
                } else { // BUGBUG: if not normalized then the varianceVector hasn't been calculated --> null pointer
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        randSigma[jjj] = dataManager.varianceVector[jjj] + ((1.0 + rand.nextDouble()) * 0.01 * dataManager.varianceVector[jjj]);
                    }
                }
                sigma[kkk] = randSigma;
            }
        }

    }
}
