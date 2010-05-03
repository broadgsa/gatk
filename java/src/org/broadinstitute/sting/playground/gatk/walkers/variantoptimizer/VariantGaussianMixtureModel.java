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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.text.XReadLines;

import Jama.*; 

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public final class VariantGaussianMixtureModel extends VariantOptimizationModel {

    protected final static Logger logger = Logger.getLogger(VariantGaussianMixtureModel.class);
    
    public final VariantDataManager dataManager;
    private final int numGaussians;
    private final int numIterations;
    private final long RANDOM_SEED = 91801305;
    private final Random rand = new Random( RANDOM_SEED );
    private final double MIN_PROB = 1E-7;
    private final double MIN_SIGMA = 1E-5;
    private final double MIN_DETERMINANT = 1E-5;

    private final double[][] mu; // The means for the clusters
    private final Matrix[] sigma; // The variances for the clusters, sigma is really sigma^2
    private final Matrix[] sigmaInverse;
    private final boolean[] deadCluster;
    private final double[] pCluster;
    private final double[] determinant;
    private final double[] clusterTITV;
    private final double[] clusterTruePositiveRate; // The true positive rate implied by the cluster's Ti/Tv ratio
    private final int minVarInCluster;
    public final boolean isUsingTiTvModel;

    private static final double INFINITE_ANNOTATION_VALUE = 6000.0;
    private static final Pattern ANNOTATION_PATTERN = Pattern.compile("^@!ANNOTATION.*");
    private static final Pattern CLUSTER_PATTERN = Pattern.compile("^@!CLUSTER.*");

    public VariantGaussianMixtureModel( final VariantDataManager _dataManager, final double _targetTITV, final int _numGaussians, final int _numIterations, final int _minVarInCluster ) {
        super( _targetTITV );
        dataManager = _dataManager;
        numGaussians = _numGaussians;
        numIterations = _numIterations;

        mu = new double[numGaussians][];
        sigma = new Matrix[numGaussians];
        determinant = new double[numGaussians];
        deadCluster = new boolean[numGaussians];
        pCluster = new double[numGaussians];
        clusterTITV = new double[numGaussians];
        clusterTruePositiveRate = new double[numGaussians];
        minVarInCluster = _minVarInCluster;
        sigmaInverse = null;
        isUsingTiTvModel = false; // this field isn't used during VariantOptimizerWalker
    }

    public VariantGaussianMixtureModel( final double _targetTITV, final String clusterFileName, final double backOffGaussianFactor ) {
        super( _targetTITV );
        final ExpandingArrayList<String> annotationLines = new ExpandingArrayList<String>();
        final ExpandingArrayList<String> clusterLines = new ExpandingArrayList<String>();

        try {
            for ( String line : new XReadLines(new File( clusterFileName )) ) {
                if( ANNOTATION_PATTERN.matcher(line).matches() ) {
                    annotationLines.add(line);
                } else if( CLUSTER_PATTERN.matcher(line).matches() ) {
                    clusterLines.add(line);
                } else {
                    throw new StingException("Malformed input file: " + clusterFileName);
                }
            }
        } catch ( FileNotFoundException e ) {
            throw new StingException("Can not find input file: " + clusterFileName);
        }

        dataManager = new VariantDataManager( annotationLines );
        // Several of the clustering parameters aren't used the second time around in ApplyVariantClusters
        numIterations = 0;
        clusterTITV = null;
        deadCluster = null;
        minVarInCluster = 0;

        // BUGBUG: move this parsing out of the constructor
        numGaussians = clusterLines.size();
        mu = new double[numGaussians][dataManager.numAnnotations];
        double sigmaVals[][][] = new double[numGaussians][dataManager.numAnnotations][dataManager.numAnnotations];
        sigma = new Matrix[numGaussians];
        sigmaInverse = new Matrix[numGaussians];
        pCluster = new double[numGaussians];
        determinant = new double[numGaussians];
        clusterTruePositiveRate = new double[numGaussians];
        boolean _isUsingTiTvModel = false;

        int kkk = 0;
        for( String line : clusterLines ) {
            final String[] vals = line.split(",");
            pCluster[kkk] = Double.parseDouble(vals[1]);
            clusterTruePositiveRate[kkk] = Double.parseDouble(vals[3]); // BUGBUG: #define these magic index numbers, very easy to make a mistake here
            if( clusterTruePositiveRate[kkk] != 1.0 ) { _isUsingTiTvModel = true; }
            for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
                mu[kkk][jjj] = Double.parseDouble(vals[4+jjj]);
                for( int ppp = 0; ppp < dataManager.numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] = Double.parseDouble(vals[4+dataManager.numAnnotations+(jjj*dataManager.numAnnotations)+ppp]) * backOffGaussianFactor; // BUGBUG: *3, suggestion by Nick to prevent GMM from over fitting and producing low likelihoods for most points
                }
            }
            
            sigma[kkk] = new Matrix(sigmaVals[kkk]);
            sigmaInverse[kkk] = sigma[kkk].inverse(); // Precompute all the inverses and determinants for use later
            determinant[kkk] = sigma[kkk].det();
            //if( determinant[kkk] < MIN_DETERMINANT ) { determinant[kkk] = MIN_DETERMINANT; }
            kkk++;
        }
        isUsingTiTvModel = _isUsingTiTvModel;

        logger.info("Found " + numGaussians + " clusters and using " + dataManager.numAnnotations + " annotations: " + dataManager.annotationKeys);
    }
    
    public final void run( final String clusterFileName ) {

        final int MAX_KNOWN_VARS = 5000000; // BUGBUG: make this a command line argument
        final int MAX_NOVEL_VARS = 5000000; // BUGBUG: make this a command line argument
        final double knownNovelMixture = 1.5; // BUGBUG: make this a command line argument

        // Create the subset of the data to cluster with
        int numNovel = 0;
        int numKnown = 0;
        for( final VariantDatum datum : dataManager.data ) {
            if( datum.isKnown ) {
                numKnown++;
            } else {
                numNovel++;
            }
        }

        final int numNovelCluster = Math.min( numNovel, MAX_NOVEL_VARS );
        final int numKnownCluster = Math.min( numKnown, MAX_KNOWN_VARS );
        final int numKnownTogether = Math.min( numKnownCluster, (int) Math.floor(knownNovelMixture * numNovelCluster) );

        final VariantDatum[] dataTogether = new VariantDatum[numNovelCluster + numKnownTogether];
        final VariantDatum[] dataKnown = new VariantDatum[numKnownCluster];

        // Create the dataTogether array, which is all the novels and 1.5x as many knowns, downsampled if there are too many
        int iii = 0;
        if( numNovelCluster == numNovel ) {
            for( final VariantDatum datum : dataManager.data ) {
                if( !datum.isKnown ) {
                    dataTogether[iii++] = datum;
                }
            }
        } else {
            logger.info("Capped at " + MAX_NOVEL_VARS + " novel variants.");
            while( iii < numNovelCluster ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( !datum.isKnown ) {
                    dataTogether[iii++] = datum;
                }
            }
        }
        if( numKnownTogether == numKnown ) {
            for( final VariantDatum datum : dataManager.data ) {
                if( datum.isKnown ) {
                    dataTogether[iii++] = datum;
                }
            }
        } else {
            while( iii < numNovelCluster + numKnownTogether ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isKnown ) {
                    dataTogether[iii++] = datum;
                }
            }
        }

        // Create the dataKnown array, which is simply all the known vars or downsampled if there are too many
        iii = 0;
        if( numKnownCluster == numKnown ) {
            for( final VariantDatum datum : dataManager.data ) {
                if( datum.isKnown ) {
                    dataKnown[iii++] = datum;
                }
            }
        } else {
            logger.info("Capped at " + MAX_KNOWN_VARS + " known variants.");
            while( iii < numKnownCluster ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isKnown ) {
                    dataKnown[iii++] = datum;
                }
            }
        }

        final boolean useTITV = true;
        logger.info("First, cluster with novels and knowns together to use ti/tv based models:");
        logger.info("Clustering with " + numNovelCluster + " novel variants and " + numKnownTogether + " known variants.");
        createClusters( dataTogether, 0, numGaussians, clusterFileName, useTITV );

        logger.info("Finally, cluster with only knowns to use ti/tv-less models:");
        logger.info("Clustering with " + numKnownCluster + " known variants.");
        createClusters( dataKnown, 0, numGaussians, clusterFileName, !useTITV );
    }

    public final void createClusters( final VariantDatum[] data, final int startCluster, final int stopCluster, final String clusterFileName, final boolean useTITV ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[][] pVarInCluster = new double[numGaussians][numVariants]; // Probability that the variant is in that cluster = simply evaluate the multivariate Gaussian

        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // ppp - loop over annotations again (full rank covariance matrix)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        // Set up the initial random Gaussians
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            deadCluster[kkk] = false;
            pCluster[kkk] = 1.0 / ((double) (stopCluster - startCluster));
            //final double[] randMu = new double[numAnnotations];
            //for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            //    randMu[jjj] = -1.5 + 3.0 * rand.nextDouble();
            //}
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            final double[][] randSigma = new double[numAnnotations][numAnnotations];
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                    randSigma[ppp][jjj] = 0.5 + 0.5 * rand.nextDouble(); // data has been normalized so sigmas are centered at 1.0
                    if(jjj != ppp) { randSigma[jjj][ppp] = 0.0; } // Sigma is a symmetric, positive-definite matrix
                }
            }
            Matrix tmp = new Matrix(randSigma);
            tmp = tmp.times(tmp.transpose());
            sigma[kkk] = tmp;
            determinant[kkk] = sigma[kkk].det();
            //if( determinant[kkk] < MIN_DETERMINANT ) { deadCluster[kkk] = true; }
        }

        // The EM loop
        for( int ttt = 0; ttt < numIterations; ttt++ ) {

            //int numValidClusters = 0;
            //for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            //    if( !deadCluster[kkk] ) { numValidClusters++; }
            //}
            //logger.info("Starting iteration " + (ttt+1) + " with " + numValidClusters + " clusters.");

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Expectation Step (calculate the probability that each data point is in each cluster)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            evaluateGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            maximizeGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Estimate each cluster's p(true) and output cluster parameters
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            outputGaussians( data, pVarInCluster, ttt+1, startCluster, stopCluster, clusterFileName, useTITV );

            logger.info("Finished iteration " + (ttt+1) );
        }
    }

    private void outputGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int iterationNumber,
                                           final int startCluster, final int stopCluster, final String clusterFileName, final boolean useTITV ) {

        if( !useTITV ) {
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                clusterTITV[kkk] = 0.0;
                clusterTruePositiveRate[kkk] = 1.0;
            }
            printClusterParamters( clusterFileName + ".WithoutTiTv." + iterationNumber );
            return;
        }

        final int numVariants = data.length;

        final double[] probTi = new double[numGaussians];
        final double[] probTv = new double[numGaussians];
        final double[] probKnown = new double[numGaussians];
        final double[] probNovel = new double[numGaussians];
        final double[] probKnownTi = new double[numGaussians];
        final double[] probKnownTv = new double[numGaussians];
        final double[] probNovelTi = new double[numGaussians];
        final double[] probNovelTv = new double[numGaussians];

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            probTi[kkk] = 0.0;
            probTv[kkk] = 0.0;
            probKnown[kkk] = 0.0;
            probNovel[kkk] = 0.0;
            probKnownTi[kkk] = 0.0;
            probKnownTv[kkk] = 0.0;
            probNovelTi[kkk] = 0.0;
            probNovelTv[kkk] = 0.0;
        }

        // Use the cluster's probabilistic Ti/Tv ratio as the indication of the cluster's true positive rate
        for( int iii = 0; iii < numVariants; iii++ ) {
            final boolean isTransition = data[iii].isTransition;
            final boolean isKnown = data[iii].isKnown;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                if( !deadCluster[kkk] ) {
                    final double prob = pVarInCluster[kkk][iii];
                    if( isKnown ) { // known
                        probKnown[kkk] += prob;
                        if( isTransition ) { // transition
                            probKnownTi[kkk] += prob;
                            probTi[kkk] += prob;
                        } else { // transversion
                            probKnownTv[kkk] += prob;
                            probTv[kkk] += prob;
                        }
                    } else { //novel
                        probNovel[kkk] += prob;
                        if( isTransition ) { // transition
                            probNovelTi[kkk] += prob;
                            probTi[kkk] += prob;
                        } else { // transversion
                            probNovelTv[kkk] += prob;
                            probTv[kkk] += prob;
                        }
                    }
                }

            }
        }

        for( int ttt = 0; ttt < 3; ttt++ ) {
            double knownAlphaFactor = 0.0;
            if( ttt == 0 ) {
                knownAlphaFactor = 0.0;
            } else if( ttt == 1 ) {
                knownAlphaFactor = 1.0;
            } else if( ttt == 2 ) {
                knownAlphaFactor = 0.5;
            }
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                if( !deadCluster[kkk] ) {
                    clusterTITV[kkk] = probTi[kkk] / probTv[kkk];
                    if( probKnown[kkk] > 500.0 && probNovel[kkk] > 500.0 ) {
                        clusterTruePositiveRate[kkk] = calcTruePositiveRateFromKnownTITV( probKnownTi[kkk] / probKnownTv[kkk], probNovelTi[kkk] / probNovelTv[kkk], clusterTITV[kkk], knownAlphaFactor );
                    } else {
                        clusterTruePositiveRate[kkk] = calcTruePositiveRateFromTITV( clusterTITV[kkk] );
                    }
                }
            }

            if( ttt == 0 ) {
                printClusterParamters( clusterFileName + ".TargetTiTv." + iterationNumber );
            } else if( ttt == 1 ) {
                printClusterParamters( clusterFileName + ".KnownTiTv." + iterationNumber );
            } else if( ttt == 2 ) {
                printClusterParamters( clusterFileName + ".BlendedTiTv." + iterationNumber );
            }
        }
    }

    private void printClusterParamters( final String clusterFileName ) {
        try {
            final PrintStream outputFile = new PrintStream( clusterFileName );
            dataManager.printClusterFileHeader( outputFile );
            final int numAnnotations = mu[0].length;
            final int numVariants = dataManager.numVariants;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( !deadCluster[kkk] ) {
                    if( pCluster[kkk] * numVariants > minVarInCluster ) {
                        final double sigmaVals[][] = sigma[kkk].getArray();
                        outputFile.print("@!CLUSTER,");
                        outputFile.print(pCluster[kkk] + ",");
                        outputFile.print(clusterTITV[kkk] + ",");
                        outputFile.print(clusterTruePositiveRate[kkk] + ",");
                        for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                            outputFile.print(mu[kkk][jjj] + ",");
                        }
                        for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                            for(int ppp = 0; ppp < numAnnotations; ppp++ ) {
                                outputFile.print(sigmaVals[jjj][ppp] + ",");
                            }
                        }
                        outputFile.println(-1);
                    }
                }
            }
            outputFile.close();
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + clusterFileName );
        }
    }

    public final double evaluateVariant( final Map<String,Object> annotationMap, final double qualityScore ) {
        final double[] pVarInCluster = new double[numGaussians];
        final double[] annotations = new double[dataManager.numAnnotations];

        for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
            double value = 0.0;
            final String annotationKey = dataManager.annotationKeys.get(jjj);
            if( annotationKey.equals("QUAL") ) {
                value = qualityScore;
            } else if( annotationKey.equals("AB") && !annotationMap.containsKey(annotationKey) ) {
                value = (0.5 - 0.005) + (0.01 * Math.random()); // HomVar calls don't have an allele balance
            } else {
                try {
                    final Object stringValue = annotationMap.get( annotationKey );
                    if( stringValue != null ) {
                        value = Double.parseDouble( stringValue.toString() );
                        if( Double.isInfinite(value) ) {
                            value = ( value > 0 ? 1.0 : -1.0 ) * INFINITE_ANNOTATION_VALUE;
                        }
                    }
                } catch( NumberFormatException e ) {
                    // do nothing, default value is 0.0
                }
            }

            annotations[jjj] = (value - dataManager.meanVector[jjj]) / dataManager.varianceVector[jjj];
        }

        evaluateGaussiansForSingleVariant( annotations, pVarInCluster );

        //if( isUsingTiTvModel ) {
            // Sum prob model
            double sum = 0.0;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                sum += pVarInCluster[kkk] * clusterTruePositiveRate[kkk];
            }
            return sum;
        /*
        } else {
            // Max prob model
            double maxProb = 0.0;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( pVarInCluster[kkk] > maxProb ) {
                    maxProb = pVarInCluster[kkk];
                }
            }
            return maxProb;
        }
        */

        // Max prob model
        /*
        double maxProb = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            if( pVarInCluster[kkk] > maxProb ) {
                maxProb = pVarInCluster[kkk];
            }
        }
        return maxProb;
        */

        // Entropy model
        /*
        double sum = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            //if( isHetCluster[kkk] == isHet ) {
                sum += pVarInCluster[kkk] * Math.log(pVarInCluster[kkk]);
            //}
        }

        double entropy = -1.0 * sum;
        double maxEntropy = -1.0 * Math.log( 1.0 / ((double) numGaussians));

        //System.out.println("H = " + entropy + ", pTrue = " + ( 1.0 - (entropy / maxEntropy) ));
        return ( 1.0 - (entropy / maxEntropy) );
        */
    }

    public final void outputOptimizationCurve( final VariantDatum[] data, final String outputPrefix, final int desiredNumVariants ) {

        final int numVariants = data.length;
        final boolean[] markedVariant = new boolean[numVariants];

        final double MAX_QUAL = 100.0;
        final double QUAL_STEP = 0.1;
        final int NUM_BINS = (int) ((MAX_QUAL / QUAL_STEP) + 1);

        final int numKnownAtCut[] = new int[NUM_BINS];
        final int numNovelAtCut[] = new int[NUM_BINS];
        final double knownTiTvAtCut[] = new double[NUM_BINS];
        final double novelTiTvAtCut[] = new double[NUM_BINS];
        final double theCut[] = new double[NUM_BINS];

        for( int iii = 0; iii < numVariants; iii++ ) {
            markedVariant[iii] = false;
        }

        PrintStream outputFile;
        try {
            outputFile = new PrintStream( outputPrefix + ".dat" );
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + outputPrefix + ".dat" );
        }

        int numKnown = 0;
        int numNovel = 0;
        int numKnownTi = 0;
        int numKnownTv = 0;
        int numNovelTi = 0;
        int numNovelTv = 0;
        boolean foundDesiredNumVariants = false;
        int jjj = 0;
        outputFile.println("pCut,numKnown,numNovel,knownTITV,novelTITV");
        for( double qCut = MAX_QUAL; qCut >= -0.001; qCut -= QUAL_STEP ) {
            for( int iii = 0; iii < numVariants; iii++ ) {
                if( !markedVariant[iii] ) {
                    if( data[iii].qual >= qCut ) {
                        markedVariant[iii] = true;
                        if( data[iii].isKnown ) { // known
                            numKnown++;
                            if( data[iii].isTransition ) { // transition
                                numKnownTi++;
                            } else { // transversion
                                numKnownTv++;
                            }
                        } else { // novel
                            numNovel++;
                            if( data[iii].isTransition ) { // transition
                                numNovelTi++;
                            } else { // transversion
                                numNovelTv++;
                            }
                        }
                    }
                }
            }
            if( desiredNumVariants != 0 && !foundDesiredNumVariants && (numKnown + numNovel) >= desiredNumVariants ) {
                logger.info( "Keeping variants with QUAL >= " + String.format("%.1f",qCut) + " results in a filtered set with: " );
                logger.info("\t" + numKnown + " known variants");
                logger.info("\t" + numNovel + " novel variants, (dbSNP rate = " + String.format("%.2f",((double) numKnown * 100.0) / ((double) numKnown + numNovel) ) + "%)");
                logger.info("\t" + String.format("%.4f known Ti/Tv ratio", ((double)numKnownTi) / ((double)numKnownTv)));
                logger.info("\t" + String.format("%.4f novel Ti/Tv ratio", ((double)numNovelTi) / ((double)numNovelTv)));
                foundDesiredNumVariants = true;
            }
            outputFile.println( qCut + "," + numKnown + "," + numNovel + "," +
                    ( numKnownTi == 0 || numKnownTv == 0 ? "NaN" : ( ((double)numKnownTi) / ((double)numKnownTv) ) ) + "," +
                    ( numNovelTi == 0 || numNovelTv == 0 ? "NaN" : ( ((double)numNovelTi) / ((double)numNovelTv) ) ));

            numKnownAtCut[jjj] = numKnown;
            numNovelAtCut[jjj] = numNovel;
            knownTiTvAtCut[jjj] = ( numKnownTi == 0 || numKnownTv == 0 ? 0.0 : ( ((double)numKnownTi) / ((double)numKnownTv) ) );
            novelTiTvAtCut[jjj] = ( numNovelTi == 0 || numNovelTv == 0 ? 0.0 : ( ((double)numNovelTi) / ((double)numNovelTv) ) );
            theCut[jjj] = qCut;
            jjj++;
        }

        // loop back through the data points looking for appropriate places to cut the data to get the target novel titv ratio
        int checkQuantile = 0;
        for( jjj = NUM_BINS-1; jjj >= 0; jjj-- ) {
            boolean foundCut = false;
            if( checkQuantile == 0 ) {
                if( novelTiTvAtCut[jjj] >= 0.9 * targetTITV ) {
                    foundCut = true;
                    checkQuantile++;
                }
            } else if( checkQuantile == 1 ) {
                if( novelTiTvAtCut[jjj] >= 0.95 * targetTITV ) {
                    foundCut = true;
                    checkQuantile++;
                }
            } else if( checkQuantile == 2 ) {
                if( novelTiTvAtCut[jjj] >= 0.98 * targetTITV ) {
                    foundCut = true;
                    checkQuantile++;
                }
            } else if( checkQuantile == 3 ) {
                if( novelTiTvAtCut[jjj] >= targetTITV ) {
                    foundCut = true;
                    checkQuantile++;
                }
            } else if( checkQuantile == 4 ) {
                break; // break out
            }

            if( foundCut ) {
                logger.info( "Keeping variants with QUAL >= " + String.format("%.1f",theCut[jjj]) + " results in a filtered set with: " );
                logger.info("\t" + numKnownAtCut[jjj] + " known variants");
                logger.info("\t" + numNovelAtCut[jjj] + " novel variants, (dbSNP rate = " +
                                    String.format("%.2f",((double) numKnownAtCut[jjj] * 100.0) / ((double) numKnownAtCut[jjj] + numNovelAtCut[jjj]) ) + "%)");
                logger.info("\t" + String.format("%.4f known Ti/Tv ratio", knownTiTvAtCut[jjj]));
                logger.info("\t" + String.format("%.4f novel Ti/Tv ratio", novelTiTvAtCut[jjj]));
            }
        }

        outputFile.close();
    }


    private void evaluateGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numAnnotations = data[0].annotations.length;
        double likelihood = 0.0;
        final double sigmaVals[][][] = new double[numGaussians][][];
        final double denom[] = new double[numGaussians];
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            if( !deadCluster[kkk] ) {
                sigmaVals[kkk] = sigma[kkk].inverse().getArray();
                denom[kkk] = Math.pow(2.0 * 3.14159, ((double)numAnnotations) / 2.0) * Math.pow(Math.abs(determinant[kkk]), 0.5);
            }
        }
        final double mult[] = new double[numAnnotations];
        for( int iii = 0; iii < data.length; iii++ ) {
            double sumProb = 0.0;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                if( !deadCluster[kkk] ) {
                    double sum = 0.0;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        mult[jjj] = 0.0;
                        for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                            mult[jjj] += (data[iii].annotations[ppp] - mu[kkk][ppp]) * sigmaVals[kkk][ppp][jjj];
                        }
                    }
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        sum += mult[jjj] * (data[iii].annotations[jjj] - mu[kkk][jjj]);
                    }

                    pVarInCluster[kkk][iii] = pCluster[kkk] * (Math.exp( -0.5 * sum ) / denom[kkk]);
                    likelihood += pVarInCluster[kkk][iii];
                    if(Double.isNaN(denom[kkk]) || determinant[kkk] < 0.5 * MIN_DETERMINANT) {
                        System.out.println("det = " + sigma[kkk].det());
                        System.out.println("denom = " + denom[kkk]);
                        System.out.println("sumExp = " + sum);
                        System.out.println("pVar = " + pVarInCluster[kkk][iii]);
                        System.out.println("=-------=");
                        throw new StingException("Numerical Instability! determinant of covariance matrix <= 0. Try running with fewer clusters and then with better behaved annotation values.");
                    }
                    if(sum < 0.0) {
                        System.out.println("det = " + sigma[kkk].det());
                        System.out.println("denom = " + denom[kkk]);
                        System.out.println("sumExp = " + sum);
                        System.out.println("pVar = " + pVarInCluster[kkk][iii]);
                        System.out.println("=-------=");
                        throw new StingException("Numerical Instability! covariance matrix no longer positive definite. Try running with fewer clusters and then with better behaved annotation values.");
                    }
                    if(pVarInCluster[kkk][iii] > 1.0) {
                        System.out.println("det = " + sigma[kkk].det());
                        System.out.println("denom = " + denom[kkk]);
                        System.out.println("sumExp = " + sum);
                        System.out.println("pVar = " + pVarInCluster[kkk][iii]);
                        System.out.println("=-------=");
                        throw new StingException("Numerical Instability! probability distribution returns > 1.0. Try running with fewer clusters and then with better behaved annotation values.");
                    }

                    if( pVarInCluster[kkk][iii] < MIN_PROB) { // Very small numbers are a very big problem
                        pVarInCluster[kkk][iii] = MIN_PROB;// + MIN_PROB * rand.nextDouble();
                    }

                    sumProb += pVarInCluster[kkk][iii];
                }
            }

            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                pVarInCluster[kkk][iii] /= sumProb;
            }

        }

        logger.info("Explained likelihood = " + String.format("%.5f",likelihood / ((double) data.length)));
    }


    private void evaluateGaussiansForSingleVariant( final double[] annotations, final double[] pVarInCluster ) {

        final int numAnnotations = annotations.length;

        double sumProb = 0.0;
        final double mult[] = new double[numAnnotations];
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            final double sigmaVals[][] = sigmaInverse[kkk].getArray();
            double sum = 0.0;
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mult[jjj] = 0.0;
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    mult[jjj] += (annotations[ppp] - mu[kkk][ppp]) * sigmaVals[ppp][jjj];
                }
            }
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                sum += mult[jjj] * (annotations[jjj] - mu[kkk][jjj]);
            }

            final double denom = Math.pow(2.0 * 3.14159, ((double)numAnnotations) / 2.0) * Math.pow(determinant[kkk], 0.5);
            pVarInCluster[kkk] =  (1.0 / ((double) numGaussians)) * (Math.exp( -0.5 * sum )) / denom;

            /*
            if( isUsingTiTvModel ) {
                //pVarInCluster[kkk] =  Math.exp( -0.5 * sum );
                if( pVarInCluster[kkk] < MIN_PROB) { // Very small numbers are a very big problem
                    pVarInCluster[kkk] = MIN_PROB;
                }
                sumProb += pVarInCluster[kkk];
            } else {
                //final double denom = Math.pow(2.0 * 3.14159, ((double)numAnnotations) / 2.0) * Math.pow(determinant[kkk], 0.5);
                //pVarInCluster[kkk] =  pCluster[kkk] * (Math.exp( -0.5 * sum )) / denom;
                //pVarInCluster[kkk] =  Math.exp( -0.5 * sum );
                // BUGBUG: should pCluster be the distribution from the GMM or a uniform distribution here?
            }
            */
        }

        //if( isUsingTiTvModel ) {
        //    for( int kkk = 0; kkk < numGaussians; kkk++ ) {
        //        pVarInCluster[kkk] /= sumProb;
        //    }
        //}
    }


    private void maximizeGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;
        final double sigmaVals[][][] = new double[numGaussians][numAnnotations][numAnnotations];

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mu[kkk][jjj] = 0.0;
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] = 0.0;
                }
            }
        }
        double sumPK = 0.0;
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            if( !deadCluster[kkk] ) {
                double sumProb = 0.0;
                for( int iii = 0; iii < numVariants; iii++ ) {
                    final double prob = pVarInCluster[kkk][iii];
                    sumProb += prob;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        mu[kkk][jjj] +=  prob * data[iii].annotations[jjj];
                    }
                }

                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    mu[kkk][jjj] /=  sumProb;
                }

                for( int iii = 0; iii < numVariants; iii++ ) {
                    final double prob = pVarInCluster[kkk][iii];
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                            sigmaVals[kkk][jjj][ppp] +=  prob * (data[iii].annotations[jjj]-mu[kkk][jjj]) * (data[iii].annotations[ppp]-mu[kkk][ppp]);
                        }
                    }
                }

                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                        if( sigmaVals[kkk][jjj][ppp] < MIN_SIGMA ) { // Very small numbers are a very big problem
                            sigmaVals[kkk][jjj][ppp] = MIN_SIGMA;// + MIN_SIGMA * rand.nextDouble();
                        }
                        sigmaVals[kkk][ppp][jjj] = sigmaVals[kkk][jjj][ppp]; // sigma must be a symmetric matrix
                    }
                }

                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        sigmaVals[kkk][jjj][ppp] /=  sumProb;
                    }
                }

                sigma[kkk] = new Matrix(sigmaVals[kkk]);
                determinant[kkk] = sigma[kkk].det();
                //if( determinant[kkk] < MIN_DETERMINANT ) { deadCluster[kkk] = true; }

                if( !deadCluster[kkk] ) {
                    pCluster[kkk] = sumProb / numVariants;
                    sumPK += pCluster[kkk];
                }
            }
        }

        // ensure pCluster sums to one, it doesn't automatically due to very small numbers getting capped
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            pCluster[kkk] /= sumPK;
        }

        /*
        // Clean up extra big or extra small clusters
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            if( pCluster[kkk] > 0.45 ) { // This is a very large cluster compared to all the others
                System.out.println("!! Found very large cluster! Busting it up into smaller clusters.");
                final int numToReplace = 3;
                final Matrix savedSigma = sigma[kkk];
                for( int rrr = 0; rrr < numToReplace; rrr++ ) {
                    // Find an example variant in the large cluster, drawn randomly
                    int randVarIndex = -1;
                    boolean foundVar = false;
                    while( !foundVar ) {
                        randVarIndex = rand.nextInt( numVariants );
                        final double probK = pVarInCluster[kkk][randVarIndex];
                        boolean inClusterK = true;
                        for( int ccc = startCluster; ccc < stopCluster; ccc++ ) {
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
                        pCluster[kkk] = 1.0 / ((double) (stopCluster-startCluster));
                    } else { // Replace the cluster with the minimum prob
                        double minProb = pCluster[startCluster];
                        int minClusterIndex = startCluster;
                        for( int ccc = startCluster; ccc < stopCluster; ccc++ ) {
                            if( pCluster[ccc] < minProb ) {
                                minProb = pCluster[ccc];
                                minClusterIndex = ccc;
                            }
                        }
                        mu[minClusterIndex] = data[randVarIndex].annotations;
                        sigma[minClusterIndex] = savedSigma;
                        //for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        //    for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        //        sigma[minClusterIndex].set(jjj, ppp, sigma[minClusterIndex].get(jjj, ppp) - 0.06 + 0.12 * rand.nextDouble());
                        //    }
                        //}
                        pCluster[minClusterIndex] = 0.5 / ((double) (stopCluster-startCluster));
                    }
                }
            }
        }
        */


        // Replace extremely small clusters with another random draw from the dataset
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            //if(determinant[kkk] < MIN_DETERMINANT ) {
            if( pCluster[kkk] < 0.0005 * (1.0 / ((double) (stopCluster-startCluster))) ||
                    determinant[kkk] < MIN_DETERMINANT ) { // This is a very small cluster compared to all the others
                logger.info("!! Found singular cluster! Initializing a new random cluster.");
                pCluster[kkk] = 0.1 / ((double) (stopCluster-startCluster)); // 0.5 /
                //final double[] randMu = new double[numAnnotations];
                //for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                //   randMu[jjj] = -1.5 + 3.0 * rand.nextDouble();
                //}
                mu[kkk] = data[rand.nextInt(numVariants)].annotations;
                final double[][] randSigma = new double[numAnnotations][numAnnotations];
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                        randSigma[ppp][jjj] = 0.50 + 0.5 * rand.nextDouble(); // data is normalized so this is centered at 1.0
                        if(jjj != ppp) { randSigma[jjj][ppp] = 0.0; } // Sigma is a symmetric, positive-definite matrix
                    }
                }
                Matrix tmp = new Matrix(randSigma);
                tmp = tmp.times(tmp.transpose());
                sigma[kkk] = tmp;
                determinant[kkk] = sigma[kkk].det();
            }
        }

        // renormalize pCluster since things might have changed due to the previous step
        sumPK = 0.0;
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            sumPK += pCluster[kkk];
        }
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            pCluster[kkk] /= sumPK;
        }
    }
}