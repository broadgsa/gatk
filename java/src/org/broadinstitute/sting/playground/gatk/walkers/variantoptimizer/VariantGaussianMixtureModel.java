package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;

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

    public final VariantDataManager dataManager;
    private final int numGaussians;
    private final int numIterations;
    private final long RANDOM_SEED = 91801305;
    private final Random rand = new Random( RANDOM_SEED );
    private final double MIN_PROB = 1E-30;
    private final double MIN_SUM_PROB = 1E-20;

    private final double[][] mu; // The means for the clusters
    private final double[][] sigma; // The variances for the clusters, sigma is really sigma^2
    private final double[] pCluster;
    //private final boolean[] isHetCluster;
    private final int[] numMaxClusterKnown;
    private final int[] numMaxClusterNovel;
    private final double[] clusterTITV;
    private final double[] clusterTruePositiveRate; // The true positive rate implied by the cluster's Ti/Tv ratio
    private final int minVarInCluster;
    private final double knownAlphaFactor;

    private static final double INFINITE_ANNOTATION_VALUE = 10000.0;
    private static final Pattern ANNOTATION_PATTERN = Pattern.compile("^@!ANNOTATION.*");
    private static final Pattern CLUSTER_PATTERN = Pattern.compile("^@!CLUSTER.*");

    public VariantGaussianMixtureModel( final VariantDataManager _dataManager, final double _targetTITV, final int _numGaussians, final int _numIterations, final int _minVarInCluster, final double _knownAlphaFactor ) {
        super( _targetTITV );
        dataManager = _dataManager;
        numGaussians = ( _numGaussians % 2 == 0 ? _numGaussians : _numGaussians + 1 );
        numIterations = _numIterations;

        mu = new double[numGaussians][];
        sigma = new double[numGaussians][];
        pCluster = new double[numGaussians];
        //isHetCluster = null;
        numMaxClusterKnown = new int[numGaussians];
        numMaxClusterNovel = new int[numGaussians];
        clusterTITV = new double[numGaussians];
        clusterTruePositiveRate = new double[numGaussians];
        minVarInCluster = _minVarInCluster;
        knownAlphaFactor = _knownAlphaFactor;
    }

    public VariantGaussianMixtureModel( final double _targetTITV, final String clusterFileName, final double backOffGaussianFactor ) {
        super( _targetTITV );
        final ExpandingArrayList<String> annotationLines = new ExpandingArrayList<String>();
        final ExpandingArrayList<String> clusterLines = new ExpandingArrayList<String>();

        try {
            for ( String line : new xReadLines(new File( clusterFileName )) ) {
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
        numIterations = 0;
        numMaxClusterKnown = null;
        numMaxClusterNovel = null;
        clusterTITV = null;
        minVarInCluster = 0;
        knownAlphaFactor = 0.0;

        //BUGBUG: move this parsing out of the constructor
        numGaussians = clusterLines.size();
        mu = new double[numGaussians][dataManager.numAnnotations];
        sigma = new double[numGaussians][dataManager.numAnnotations];
        pCluster = new double[numGaussians];
        //isHetCluster = new boolean[numGaussians];
        clusterTruePositiveRate = new double[numGaussians];

        int kkk = 0;
        for( String line : clusterLines ) {
            final String[] vals = line.split(",");
            //isHetCluster[kkk] = Integer.parseInt(vals[1]) == 1;
            pCluster[kkk] = Double.parseDouble(vals[2]);
            clusterTruePositiveRate[kkk] = Double.parseDouble(vals[6]); //BUGBUG: #define these magic index numbers, very easy to make a mistake here
            for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
                mu[kkk][jjj] = Double.parseDouble(vals[7+jjj]);
                sigma[kkk][jjj] = Double.parseDouble(vals[7+dataManager.numAnnotations+jjj]) * backOffGaussianFactor; //BUGBUG: *3, suggestion by Nick to prevent GMM from over fitting and producing low likelihoods for most points
            }
            kkk++;
        }

        System.out.println("Found " + numGaussians + " clusters and using " + dataManager.numAnnotations + " annotations: " + dataManager.annotationKeys);
    }
    
    public final void run( final String clusterFileName ) {

        final int MAX_VARS = 1000000; //BUGBUG: make this a command line argument

        // Create the subset of the data to cluster with
        int numNovel = 0;
        int numKnown = 0;
        int numHet = 0;
        int numHom = 0;
        for( final VariantDatum datum : dataManager.data ) {
            if( datum.isKnown ) {
                numKnown++;
            } else {
                numNovel++;
            }
            if( datum.isHet ) {
                numHet++;
            } else {
                numHom++;
            }
        }

        // This block of code is used to cluster with novels + 1.5x knowns mixed together

        VariantDatum[] data;

        // Grab a set of data that is all of the novel variants plus 1.5x as many known variants drawn at random
        // If there are almost as many novels as known, simply use all the variants
        // BUGBUG: allow downsampling and arbitrary mixtures of knowns and novels
        final int numSubset = (int)Math.floor(numNovel*2.5);
        if( numSubset * 1.3 < dataManager.numVariants ) {
            data = new VariantDatum[numSubset];
            int iii = 0;
            for( final VariantDatum datum : dataManager.data ) {
                if( !datum.isKnown ) {
                    data[iii++] = datum;
                }
            }
            while( iii < numSubset ) { // grab an equal number of known variants at random
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isKnown ) {
                    data[iii++] = datum;
                }
            }
        } else {
            data = dataManager.data;
        }

        System.out.println("Clustering with " + numNovel + " novel variants and " + (data.length - numNovel) + " known variants...");
        if( data.length == dataManager.numVariants ) { System.out.println(" (used all variants since 2.5*numNovel is so large compared to the full set) "); }
        createClusters( data, 0, numGaussians ); // Using a subset of the data
        System.out.println("Outputting cluster parameters...");
        printClusters( clusterFileName );





        // This block of code is to cluster knowns and novels separately
        /*
        final VariantDatum[] dataNovel = new VariantDatum[Math.min(numNovel,MAX_VARS)];
        final VariantDatum[] dataKnown = new VariantDatum[Math.min(numKnown,MAX_VARS)];

        //BUGBUG: This is ugly
        int jjj = 0;
        if(numNovel <= MAX_VARS) {
            for( final VariantDatum datum : dataManager.data ) {
                if( !datum.isKnown ) {
                    dataNovel[jjj++] = datum;
                }
            }
        } else {
            System.out.println("Capped at " + MAX_VARS + " novel variants.");
            while( jjj < MAX_VARS ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( !datum.isKnown ) {
                    dataNovel[jjj++] = datum;
                }
            }
        }

        int iii = 0;
        if(numKnown <= MAX_VARS) {
            for( final VariantDatum datum : dataManager.data ) {
                if( datum.isKnown ) {
                    dataKnown[iii++] = datum;
                }
            }
        } else {
            System.out.println("Capped at " + MAX_VARS + " known variants.");
            while( iii < MAX_VARS ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isKnown ) {
                    dataKnown[iii++] = datum;
                }
            }
        }


        System.out.println("Clustering with " + Math.min(numNovel,MAX_VARS) + " novel variants.");
        createClusters( dataNovel, 0, numGaussians / 2 );
        System.out.println("Clustering with " + Math.min(numKnown,MAX_VARS) + " known variants.");
        createClusters( dataKnown, numGaussians / 2, numGaussians );
        System.out.println("Outputting cluster parameters...");
        printClusters( clusterFileName );

        */

        /*
        // This block of code is to cluster het and hom calls separately, but mixing together knowns and novels
        final VariantDatum[] dataHet = new VariantDatum[Math.min(numHet,MAX_VARS)];
        final VariantDatum[] dataHom = new VariantDatum[Math.min(numHom,MAX_VARS)];

        //BUGBUG: This is ugly
        int jjj = 0;
        if(numHet <= MAX_VARS) {
            for( final VariantDatum datum : dataManager.data ) {
                if( datum.isHet ) {
                    dataHet[jjj++] = datum;
                }
            }
        } else {
            System.out.println("Found " + numHet + " het variants but capped at clustering with " + MAX_VARS + ".");
            while( jjj < MAX_VARS ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( datum.isHet ) {
                    dataHet[jjj++] = datum;
                }
            }
        }

        int iii = 0;
        if(numHom <= MAX_VARS) {
            for( final VariantDatum datum : dataManager.data ) {
                if( !datum.isHet ) {
                    dataHom[iii++] = datum;
                }
            }
        } else {
            System.out.println("Found " + numHom + " hom variants but capped at clustering with " + MAX_VARS + ".");
            while( iii < MAX_VARS ) {
                final VariantDatum datum = dataManager.data[rand.nextInt(dataManager.numVariants)];
                if( !datum.isHet ) {
                    dataHom[iii++] = datum;
                }
            }
        }

        System.out.println("Clustering with " + Math.min(numHet,MAX_VARS) + " het variants.");
        createClusters( dataHet, 0, numGaussians / 2 );
        System.out.println("Clustering with " + Math.min(numHom,MAX_VARS) + " hom variants.");
        createClusters( dataHom, numGaussians / 2, numGaussians );
        System.out.println("Outputting cluster parameters...");
        printClusters( clusterFileName );
*/
    }


/*
    public final void createClusters( final VariantDatum[] data, int startCluster, int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[][] pVarInCluster = new double[numGaussians][numVariants];
        final double[] probTi = new double[numGaussians];
        final double[] probTv = new double[numGaussians];
        final double[] probKnown = new double[numGaussians];
        final double[] probNovel = new double[numGaussians];
        final double[] probKnownTi = new double[numGaussians];
        final double[] probKnownTv = new double[numGaussians];
        final double[] probNovelTi = new double[numGaussians];
        final double[] probNovelTv = new double[numGaussians];

        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        // Set up the initial random Gaussians
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            numMaxClusterKnown[kkk] = 0;
            numMaxClusterNovel[kkk] = 0;
            pCluster[kkk] = 1.0 / ((double) (stopCluster - startCluster));
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            final double[] randSigma = new double[numAnnotations];
            if( dataManager.isNormalized ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = 0.7 + 0.4 * rand.nextDouble();
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
            evaluateGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            maximizeGaussians( data, pVarInCluster, startCluster, stopCluster );

            System.out.println("Finished iteration " + (ttt+1) );
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Evaluate the clusters using titv as an estimate of the true positive rate
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        evaluateGaussians( data, pVarInCluster, startCluster, stopCluster ); // One final evaluation because the Gaussians moved in the last maximization step

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            probTi[kkk] = 0.0;
            probTv[kkk] = 0.0;
            probKnown[kkk] = 0.0;
            probNovel[kkk] = 0.0;
        }
        for( int iii = 0; iii < numVariants; iii++ ) {
            final boolean isTransition = data[iii].isTransition;
            final boolean isKnown = data[iii].isKnown;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
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

            double maxProb = pVarInCluster[startCluster][iii];
            int maxCluster = startCluster;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                if( pVarInCluster[kkk][iii] > maxProb ) {
                    maxProb = pVarInCluster[kkk][iii];
                    maxCluster = kkk;
                }
            }
            if( isKnown ) {
                numMaxClusterKnown[maxCluster]++;
            } else {
                numMaxClusterNovel[maxCluster]++;
            }
        }

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            clusterTITV[kkk] = probTi[kkk] / probTv[kkk];
            if( probKnown[kkk] > 600.0 ) { // BUGBUG: make this a command line argument, parameterize performance based on this important argument
                clusterTruePositiveRate[kkk] = calcTruePositiveRateFromKnownTITV( probKnownTi[kkk] / probKnownTv[kkk], probNovelTi[kkk] / probNovelTv[kkk], clusterTITV[kkk], knownAlphaFactor );
            } else {
                clusterTruePositiveRate[kkk] = calcTruePositiveRateFromTITV( clusterTITV[kkk] );
            }
        }
    }
*/


    // This cluster method doesn't make use of the differences between known and novel Ti/Tv ratios


    public final void createClusters( final VariantDatum[] data, final int startCluster, final int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[][] pVarInCluster = new double[numGaussians][numVariants]; // Probability that the variant is in that cluster = simply evaluate the multivariate Gaussian
        final double[] probTi = new double[numGaussians];
        final double[] probTv = new double[numGaussians];

        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        // Set up the initial random Gaussians
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            numMaxClusterKnown[kkk] = 0;
            numMaxClusterNovel[kkk] = 0;
            pCluster[kkk] = 1.0 / ((double) (stopCluster - startCluster));
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            final double[] randSigma = new double[numAnnotations];
            if( dataManager.isNormalized ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = 0.7 + 0.4 * rand.nextDouble();
                }
            } else { // BUGBUG: if not normalized then the varianceVector hasn't been calculated --> null pointer
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    randSigma[jjj] = dataManager.varianceVector[jjj] + ((1.0 + rand.nextDouble()) * 0.01 * dataManager.varianceVector[jjj]);
                }
            }
            sigma[kkk] = randSigma;
        }

        // The EM loop
        for( int ttt = 0; ttt < numIterations; ttt++ ) {

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Expectation Step (calculate the probability that each data point is in each cluster)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            evaluateGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            maximizeGaussians( data, pVarInCluster, startCluster, stopCluster );

            System.out.println("Finished iteration " + (ttt+1) );
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Evaluate the clusters using titv as an estimate of the true positive rate
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        evaluateGaussians( data, pVarInCluster, startCluster, stopCluster ); // One final evaluation because the Gaussians moved in the last maximization step

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            probTi[kkk] = 0.0;
            probTv[kkk] = 0.0;
        }
        // Use the cluster's probabilistic Ti/Tv ratio as the indication of the cluster's true positive rate
        for( int iii = 0; iii < numVariants; iii++ ) {
            if( data[iii].isTransition ) { // transition
                for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                    probTi[kkk] += pVarInCluster[kkk][iii];
                }
            } else { // transversion
                for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                    probTv[kkk] += pVarInCluster[kkk][iii];
                }
            }

            // Calculate which cluster has the maximum probability for this variant for use as a metric of how well clustered the data is 
            double maxProb = pVarInCluster[startCluster][iii];
            int maxCluster = startCluster;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
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
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            clusterTITV[kkk] = probTi[kkk] / probTv[kkk];
            clusterTruePositiveRate[kkk] = calcTruePositiveRateFromTITV( clusterTITV[kkk] );
        }
    }



    private void printClusters( final String clusterFileName ) {
        try {
            final PrintStream outputFile = new PrintStream( clusterFileName );
            dataManager.printClusterFileHeader( outputFile );
            final int numAnnotations = mu[0].length;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( numMaxClusterKnown[kkk] + numMaxClusterNovel[kkk] >= minVarInCluster ) {
                    outputFile.print("@!CLUSTER,");
                    outputFile.print( (kkk < numGaussians / 2 ? 1 : 0) + "," ); // is het cluster?
                    outputFile.print(pCluster[kkk] + ",");
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
                }
            }
            outputFile.close();
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + clusterFileName );
        }
    }

    public final double evaluateVariant( final Map<String,String> annotationMap, final double qualityScore, final boolean isHet ) {
        final double[] pVarInCluster = new double[numGaussians];
        final double[] annotations = new double[dataManager.numAnnotations];

        for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
            double value = 0.0;
            if( dataManager.annotationKeys.get(jjj).equals("QUAL") ) {
                value = qualityScore;
            } else {
                try {
                    final String stringValue = annotationMap.get( dataManager.annotationKeys.get(jjj) );
                    if( stringValue != null ) {
                        value = Double.parseDouble( stringValue );
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

        evaluateGaussiansForSingleVariant( annotations, pVarInCluster, isHet );

        double sum = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            //if( isHetCluster[kkk] == isHet ) {
                sum += pVarInCluster[kkk] * clusterTruePositiveRate[kkk];
            //}
        }

        return sum;
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
        for( double qCut = MAX_QUAL; qCut >= 0.0; qCut -= QUAL_STEP ) {
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
                System.out.println( "Keeping variants with QUAL >= " + String.format("%.1f",qCut) + " results in a filtered set with: " );
                System.out.println("\t" + numKnown + " known variants");
                System.out.println("\t" + numNovel + " novel variants, (dbSNP rate = " + String.format("%.2f",((double) numKnown * 100.0) / ((double) numKnown + numNovel) ) + "%)");
                System.out.println("\t" + String.format("%.4f known Ti/Tv ratio", ((double)numKnownTi) / ((double)numKnownTv)));
                System.out.println("\t" + String.format("%.4f novel Ti/Tv ratio", ((double)numNovelTi) / ((double)numNovelTv)));
                System.out.println();
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
                System.out.println( "Keeping variants with QUAL >= " + String.format("%.1f",theCut[jjj]) + " results in a filtered set with: " );
                System.out.println("\t" + numKnownAtCut[jjj] + " known variants");
                System.out.println("\t" + numNovelAtCut[jjj] + " novel variants, (dbSNP rate = " +
                                    String.format("%.2f",((double) numKnownAtCut[jjj] * 100.0) / ((double) numKnownAtCut[jjj] + numNovelAtCut[jjj]) ) + "%)");
                System.out.println("\t" + String.format("%.4f known Ti/Tv ratio", knownTiTvAtCut[jjj]));
                System.out.println("\t" + String.format("%.4f novel Ti/Tv ratio", novelTiTvAtCut[jjj]));
                System.out.println();
            }
        }

        outputFile.close();
    }


    private void evaluateGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numAnnotations = data[0].annotations.length;

        for( int iii = 0; iii < data.length; iii++ ) {
            double sumProb = 0.0;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
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
                for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                    pVarInCluster[kkk][iii] /= sumProb;
                }
            }
        }
    }


    private void evaluateGaussiansForSingleVariant( final double[] annotations, final double[] pVarInCluster, final boolean isHet ) {

        final int numAnnotations = annotations.length;

        double sumProb = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            //if( isHetCluster[kkk] == isHet ) {
                double sum = 0.0;
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    sum += ( (annotations[jjj] - mu[kkk][jjj]) * (annotations[jjj] - mu[kkk][jjj]) )
                            / sigma[kkk][jjj];
                }

                // BUGBUG: reverting to old version that didn't have pCluster[kkk]* here, this meant that the overfitting parameters changed meanings
                //pVarInCluster[kkk] = pCluster[kkk] * Math.exp( -0.5 * sum );
                pVarInCluster[kkk] = Math.exp( -0.5 * sum );

                if( pVarInCluster[kkk] < MIN_PROB) { // Very small numbers are a very big problem
                    pVarInCluster[kkk] = MIN_PROB;
                }

                sumProb += pVarInCluster[kkk];
            //}
        }

        if( sumProb > MIN_SUM_PROB ) { // Very small numbers are a very big problem
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                //if( isHetCluster[kkk] == isHet ) {
                    pVarInCluster[kkk] /= sumProb;
                //}
            }
        }
    }


    private void maximizeGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mu[kkk][jjj] = 0.0;
                sigma[kkk][jjj] = 0.0;
            }
        }
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
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
            } //BUGBUG: clean up dead clusters to speed up computation

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
        //BUGBUG: Is this a good idea?
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            if( pCluster[kkk] > 0.45 ) { // This is a very large cluster compared to all the others
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
                        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                            sigma[minClusterIndex][jjj] += -0.06 + 0.12 * rand.nextDouble();
                            if( sigma[minClusterIndex][jjj] < MIN_SUM_PROB ) {
                                sigma[minClusterIndex][jjj] = MIN_SUM_PROB;
                            }
                        }
                        pCluster[minClusterIndex] = 1.0 / ((double) (stopCluster-startCluster));
                    }
                }
            }
        }
        

        // Replace small clusters with another random draw from the dataset
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            if( pCluster[kkk] < 0.07 * (1.0 / ((double) (stopCluster-startCluster))) ) { // This is a very small cluster compared to all the others
                pCluster[kkk] = 1.0 / ((double) (stopCluster-startCluster));
                mu[kkk] = data[rand.nextInt(numVariants)].annotations;
                final double[] randSigma = new double[numAnnotations];
                if( dataManager.isNormalized ) {
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        randSigma[jjj] = 0.7 + 0.4 * rand.nextDouble(); // BUGBUG: Explore a wider range of possible sigma values since we are tossing out clusters anyway?
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