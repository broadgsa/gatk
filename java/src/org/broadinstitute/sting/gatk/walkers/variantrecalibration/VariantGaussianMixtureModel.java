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

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.text.XReadLines;

import Jama.*; 

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
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
    private final int maxIterations;
    private final static long RANDOM_SEED = 91801305;
    private final static Random rand = new Random( RANDOM_SEED );
    private final double MIN_SIGMA = 1E-10;
    private final double MIN_DETERMINANT = 1E-3;
    private final double MIN_PROB_CONVERGENCE = 1E-5;
    private final boolean FORCE_INDEPENDENT_ANNOTATIONS;

    private final double[][] mu; // The means for each cluster
    private final Matrix[] sigma; // The covariance matrix for each cluster
    private final Matrix[] sigmaInverse;
    private double[] pClusterLog10;
    private final double[] determinant;
    private final double[] alleleCountFactorArray;
    private final int minVarInCluster;
    private final double stdThreshold;
    private final double qualThreshold;

    private static final Pattern ANNOTATION_PATTERN = Pattern.compile("^@!ANNOTATION.*");
    private static final Pattern ALLELECOUNT_PATTERN = Pattern.compile("^@!ALLELECOUNT.*");
    private static final Pattern CLUSTER_PATTERN = Pattern.compile("^@!CLUSTER.*");

    public VariantGaussianMixtureModel( final VariantDataManager _dataManager, final int _numGaussians, final int _maxIterations, final int _minVarInCluster,
                                        final int maxAC, final boolean _forceIndependent, final double _stdThreshold, final double _qualThreshold ) {
        dataManager = _dataManager;
        numGaussians = _numGaussians;
        maxIterations = _maxIterations;

        mu = new double[numGaussians][];
        sigma = new Matrix[numGaussians];
        determinant = new double[numGaussians];
        pClusterLog10 = new double[numGaussians];
        alleleCountFactorArray = new double[maxAC + 1];
        minVarInCluster = _minVarInCluster;
        stdThreshold = _stdThreshold;
        qualThreshold = _qualThreshold;
        FORCE_INDEPENDENT_ANNOTATIONS = _forceIndependent;
        sigmaInverse = null; // This field isn't used during VariantOptimizer pass
    }

    public VariantGaussianMixtureModel( final double _targetTITV, final String clusterFileName, final double backOffGaussianFactor ) {
        super( _targetTITV );
        final ExpandingArrayList<String> annotationLines = new ExpandingArrayList<String>();
        final ExpandingArrayList<String> alleleCountLines = new ExpandingArrayList<String>();
        final ExpandingArrayList<String> clusterLines = new ExpandingArrayList<String>();

        try {
            for ( final String line : new XReadLines(new File( clusterFileName )) ) {
                if( ANNOTATION_PATTERN.matcher(line).matches() ) {
                    annotationLines.add(line);
                } else if( ALLELECOUNT_PATTERN.matcher(line).matches() ) {
                    alleleCountLines.add(line);
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
        maxIterations = 0;
        minVarInCluster = 0;
        stdThreshold = 0.0;
        qualThreshold = 0.0;
        FORCE_INDEPENDENT_ANNOTATIONS = false;

        // BUGBUG: move this parsing out of the constructor
        numGaussians = clusterLines.size();
        mu = new double[numGaussians][dataManager.numAnnotations];
        final double sigmaVals[][][] = new double[numGaussians][dataManager.numAnnotations][dataManager.numAnnotations];
        sigma = new Matrix[numGaussians];
        sigmaInverse = new Matrix[numGaussians];
        pClusterLog10 = new double[numGaussians];
        determinant = new double[numGaussians];

        alleleCountFactorArray = new double[alleleCountLines.size() + 1];
        for( final String line : alleleCountLines ) {
            final String[] vals = line.split(",");
            alleleCountFactorArray[Integer.parseInt(vals[1])] = Double.parseDouble(vals[2]);
        }

        int kkk = 0;
        for( final String line : clusterLines ) {
            final String[] vals = line.split(",");
            pClusterLog10[kkk] = Math.log10( Double.parseDouble(vals[1]) ); // BUGBUG: #define these magic index numbers, very easy to make a mistake here
            for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
                mu[kkk][jjj] = Double.parseDouble(vals[2+jjj]);
                for( int ppp = 0; ppp < dataManager.numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] = Double.parseDouble(vals[2+dataManager.numAnnotations+(jjj*dataManager.numAnnotations)+ppp]) * backOffGaussianFactor;
                }
            }
            
            sigma[kkk] = new Matrix(sigmaVals[kkk]);
            sigmaInverse[kkk] = sigma[kkk].inverse(); // Precompute all the inverses and determinants for use later
            determinant[kkk] = sigma[kkk].det();
            kkk++;
        }

        logger.info("Found " + numGaussians + " clusters using " + dataManager.numAnnotations + " annotations: " + dataManager.annotationKeys);
    }
    
    public final void run( final String clusterFileName ) {

        // Initialize the Allele Count prior
        generateAlleleCountPrior();

        int numValid = 0;
        int numOutlier = 0;
        int numBadQual = 0;
        int numZeroWeight = 0;

        // Only cluster with a good set of knowns. Filter based on being too many std's away from the mean annotation value or having low quality score
        for( final VariantDatum datum : dataManager.data ) {
            boolean goodVar = true;
            if(!(datum.weight > 0.0)) {
                goodVar = false;
                numZeroWeight++;
            }
            if(goodVar) {
                for( final double val : datum.annotations ) {
                    if( Math.abs(val) > stdThreshold ) {
                        goodVar = false;
                        numOutlier++;
                        break;
                    }
                }
            }
            if(goodVar) {
                if( datum.qual < qualThreshold ) {
                    goodVar = false;
                    numBadQual++;
                }
            }
            if(goodVar) { numValid++; }
        }

        final VariantDatum data[] = new VariantDatum[numValid];
        int iii = 0;
        for( final VariantDatum datum : dataManager.data ) {
            boolean goodVar = true;
            if(!(datum.weight > 0.0)) {
                goodVar = false;
            }
            if(goodVar) {
                for( final double val : datum.annotations ) {
                    if( Math.abs(val) > stdThreshold ) {
                        goodVar = false;
                        break;
                    }
                }
            }
            if(goodVar) {
                if( datum.qual < qualThreshold ) {
                    goodVar = false;
                }
            }
            if(goodVar) { data[iii++] = datum; }
        }

        logger.info("Clustering with " + data.length + " valid variants.");
        logger.info("  " + numZeroWeight + " variants were removed from clustering due to having zero clustering weight.");
        logger.info("  " + numOutlier + " variants were removed due to having annotations that were more than " + stdThreshold + " standard deviations away from the mean annotation value.");
        logger.info("  " + numBadQual + " variants were removed because raw QUAL value was less than threshold (" + qualThreshold + ").");
        createClusters( data, 0, numGaussians, clusterFileName );

        // Simply cluster with all the variants. The knowns have been given more weight than the novels
        //logger.info("Clustering with " + dataManager.data.length + " variants.");
        //createClusters( dataManager.data, 0, numGaussians, clusterFileName );
    }

    private void generateAlleleCountPrior() {

        final double[] acExpectation = new double[alleleCountFactorArray.length];
        final double[] acActual = new double[alleleCountFactorArray.length];
        final int[] alleleCount = new int[alleleCountFactorArray.length];

        double sumExpectation = 0.0;
        for( int iii = 1; iii < alleleCountFactorArray.length; iii++ ) {
            acExpectation[iii] = 1.0 / ((double) iii);
            sumExpectation += acExpectation[iii];
        }
        for( int iii = 1; iii < alleleCountFactorArray.length; iii++ ) {
            acExpectation[iii] /= sumExpectation; // Turn acExpectation into a probability distribution
            alleleCount[iii] = 1; // Start off with one count to smooth the estimate
        }
        for( final VariantDatum datum : dataManager.data ) {
            alleleCount[datum.alleleCount]++;
        }
        for( int iii = 1; iii < alleleCountFactorArray.length; iii++ ) {
            acActual[iii] = ((double)alleleCount[iii]) / ((double) (dataManager.data.length+(alleleCountFactorArray.length-1))); // Turn acActual into a probability distribution
        }
        for( int iii = 1; iii < alleleCountFactorArray.length; iii++ ) {
            alleleCountFactorArray[iii] = acExpectation[iii] / acActual[iii]; // Prior is (expected / observed)
        }
    }

    public final double getAlleleCountPrior( final int alleleCount ) {
        return alleleCountFactorArray[alleleCount];
    }

    public final void createClusters( final VariantDatum[] data, final int startCluster, final int stopCluster, final String clusterFileName ) {

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
            pClusterLog10[kkk] = Math.log10(1.0 / ((double) (stopCluster - startCluster)));
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            final double[][] randSigma = new double[numAnnotations][numAnnotations];
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                    randSigma[ppp][jjj] = 0.55 + 1.25 * rand.nextDouble();
                    if(rand.nextBoolean()) {
                        randSigma[ppp][jjj] *= -1.0;
                    }
                    if(jjj != ppp) { randSigma[jjj][ppp] = 0.0; } // Sigma is a symmetric, positive-definite matrix created by taking a lower diagonal matrix and multiplying by its transpose
                }
            }
            if( FORCE_INDEPENDENT_ANNOTATIONS ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        if(jjj!=ppp) {
                            randSigma[jjj][ppp] = 0.0;
                        }
                    }
                }

            }
            Matrix tmp = new Matrix(randSigma);
            tmp = tmp.times(tmp.transpose());
            sigma[kkk] = tmp;
            determinant[kkk] = sigma[kkk].det();
        }

        // The EM loop
        double previousLikelihood = -1E20;
        double currentLikelihood;
        int ttt = 1;
        while( ttt < maxIterations ) {

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Expectation Step (calculate the probability that each data point is in each cluster)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            currentLikelihood = evaluateGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Maximization Step (move the clusters to maximize the sum probability of each data point)
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            maximizeGaussians( data, pVarInCluster, startCluster, stopCluster );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Output cluster parameters at each iteration
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            printClusterParameters( clusterFileName + "." + ttt );

            logger.info("Finished iteration " + ttt );
            ttt++;
            if( currentLikelihood - previousLikelihood < MIN_PROB_CONVERGENCE) {
                logger.info("Convergence!");
                break;
            }
            previousLikelihood = currentLikelihood;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Output the final cluster parameters
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printClusterParameters( clusterFileName );
    }

    private void printClusterParameters( final String clusterFileName ) {
        try {
            final PrintStream outputFile = new PrintStream( clusterFileName );
            dataManager.printClusterFileHeader( outputFile );
            for( int iii = 1; iii < alleleCountFactorArray.length; iii++ ) {
                outputFile.print("@!ALLELECOUNT,");
                outputFile.println(iii + "," + alleleCountFactorArray[iii]);
            }

            final int numAnnotations = mu[0].length;
            final int numVariants = dataManager.numVariants;
            for( int kkk = 0; kkk < numGaussians; kkk++ ) {
                if( Math.pow(10.0, pClusterLog10[kkk]) * numVariants > minVarInCluster ) {
                    final double sigmaVals[][] = sigma[kkk].getArray();
                    outputFile.print("@!CLUSTER,");
                    outputFile.print(Math.pow(10.0, pClusterLog10[kkk]) + ",");
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
            outputFile.close();
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + clusterFileName );
        }
    }

    public static double decodeAnnotation( final String annotationKey, final VariantContext vc ) {
        double value;
        //if( annotationKey.equals("AB") && !vc.getAttributes().containsKey(annotationKey) ) {
        //    value = (0.5 - 0.005) + (0.01 * rand.nextDouble()); // HomVar calls don't have an allele balance
        //}
        if( annotationKey.equals("QUAL") ) {
            value = vc.getPhredScaledQual();
        } else {
            try {
                value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            } catch( Exception e ) {
                throw new StingException("No double value detected for annotation = " + annotationKey +
                        " in variant at " + vc.getLocation() + ", reported annotation value = " + vc.getAttribute( annotationKey ) ); 
            }
        }
        return value;
    }

    public final double evaluateVariantLog10( final VariantContext vc ) {
        final double[] pVarInCluster = new double[numGaussians];
        final double[] annotations = new double[dataManager.numAnnotations];

        for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
            final double value = decodeAnnotation( dataManager.annotationKeys.get(jjj), vc );
            annotations[jjj] = (value - dataManager.meanVector[jjj]) / dataManager.varianceVector[jjj];
        }

        evaluateGaussiansForSingleVariant( annotations, pVarInCluster );

        double sum = 0.0;
        for( int kkk = 0; kkk < numGaussians; kkk++ ) {
            sum += pVarInCluster[kkk]; // * clusterTruePositiveRate[kkk];
        }

        double sumLog10 = Math.log10(sum);
        if( Double.isInfinite(sumLog10) || Double.isNaN(sumLog10) ) {
            //logger.warn("pTrueLog10 = -Infinity, capped at -20");
            sumLog10 = -20;
        }
        return sumLog10;
    }

   public final void outputClusterReports( final String outputPrefix ) {
        final double STD_STEP = 0.2;
        final double MAX_STD = 4.0;
        final double MIN_STD = -4.0;
        final int NUM_BINS = (int)Math.floor((Math.abs(MIN_STD) + Math.abs(MAX_STD)) / STD_STEP);
        final int numAnnotations = dataManager.numAnnotations;
        int totalCountsKnown = 0;
        int totalCountsNovel = 0;

        final int counts[][][] = new int[numAnnotations][NUM_BINS][2];
        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                counts[jjj][iii][0] = 0;
                counts[jjj][iii][1] = 0;
            }
        }

        for( final VariantDatum datum : dataManager.data ) {
            final int isKnown = ( datum.isKnown ? 1 : 0 );
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                int histBin = (int)Math.round((datum.annotations[jjj]-MIN_STD) * (1.0 / STD_STEP));
                if(histBin < 0) { histBin = 0; }
                if(histBin > NUM_BINS-1) { histBin = NUM_BINS-1; }
                if(histBin >= 0 && histBin <= NUM_BINS-1) {
                    counts[jjj][histBin][isKnown]++;
                }
            }
            if( isKnown == 1 ) { totalCountsKnown++; }
            else { totalCountsNovel++; }
        }

        int annIndex = 0;
        for( final String annotation : dataManager.annotationKeys ) {
            PrintStream outputFile;
            try {
                outputFile = new PrintStream( outputPrefix + "." + annotation + ".dat" );
            } catch (Exception e) {
                throw new StingException( "Unable to create output file: " + outputPrefix + ".dat" );
            }

            outputFile.println("annotationValue,knownDist,novelDist");

            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                final double annotationValue = (((double)iii * STD_STEP)+MIN_STD) * dataManager.varianceVector[annIndex] + dataManager.meanVector[annIndex];
                outputFile.println( annotationValue + "," + ( ((double)counts[annIndex][iii][1])/((double)totalCountsKnown) ) +
                                                      "," + ( ((double)counts[annIndex][iii][0])/((double)totalCountsNovel) ));
            }

            annIndex++;
        }

       // BUGBUG: next output the actual cluster on top by integrating out every other annotation
    }

    public final void outputOptimizationCurve( final VariantDatum[] data, final String outputPrefix,
                                               final int desiredNumVariants, final Double[] FDRtranches, final double QUAL_STEP ) {

        final int numVariants = data.length;
        final boolean[] markedVariant = new boolean[numVariants];

        final double MAX_QUAL = 100.0;
        final int NUM_BINS = (int) ((MAX_QUAL / QUAL_STEP) + 1);

        final int numKnownAtCut[] = new int[NUM_BINS];
        final int numNovelAtCut[] = new int[NUM_BINS];
        final double knownTiTvAtCut[] = new double[NUM_BINS];
        final double novelTiTvAtCut[] = new double[NUM_BINS];
        final double theCut[] = new double[NUM_BINS];

        final double fdrCutAsTiTv[] = new double[FDRtranches.length];
        for( int iii = 0; iii < FDRtranches.length; iii++ ) {
            fdrCutAsTiTv[iii] = (1.0 - FDRtranches[iii] / 100.0) * (targetTITV - 0.5) + 0.5;
        }

        for( int iii = 0; iii < numVariants; iii++ ) {
            markedVariant[iii] = false;
        }

        PrintStream outputFile;
        try {
            outputFile = new PrintStream( outputPrefix + ".dat" );
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + outputPrefix + ".dat" );
        }
        PrintStream tranchesOutputFile;
        try {
            tranchesOutputFile = new PrintStream( outputPrefix + ".dat.tranches" );
            tranchesOutputFile.println("FDRtranche,novelTITV,pCut,numNovel,filterName");
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + outputPrefix + ".dat.tranches" );
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
        int tranche = FDRtranches.length - 1;
        for( jjj = NUM_BINS-1; jjj >= 0; jjj-- ) {

            if( tranche >= 0 && novelTiTvAtCut[jjj] >= fdrCutAsTiTv[tranche] ) {
                tranchesOutputFile.println(String.format("%.2f,%.2f,%.2f,%d,FDRtranche%.2fto%.2f",
                        FDRtranches[tranche],novelTiTvAtCut[jjj],theCut[jjj],numNovelAtCut[jjj],
                        (tranche == 0 ? 0.0 : FDRtranches[tranche-1]) ,FDRtranches[tranche]));
                tranche--;
            }

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
                logger.info("\t" + String.format("--> with an implied novel FDR of %.2f percent", Math.abs(100.0 * (1.0-((novelTiTvAtCut[jjj] - 0.5) / (targetTITV - 0.5))))));
            }
        }

        outputFile.close();
        tranchesOutputFile.close();
    }


    private double evaluateGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numAnnotations = data[0].annotations.length;
        double likelihood = 0.0;
        final double sigmaVals[][][] = new double[numGaussians][][];
        final double denomLog10[] = new double[numGaussians];
        final double pVarInClusterLog10[] = new double[numGaussians];
        double pVarInClusterReals[];

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            sigmaVals[kkk] = sigma[kkk].inverse().getArray();
            denomLog10[kkk] = Math.log10(Math.pow(2.0 * Math.PI, ((double)numAnnotations) / 2.0)) + Math.log10(Math.pow(determinant[kkk], 0.5));
            if( Double.isInfinite(denomLog10[kkk]) ) {
                throw new StingException("Numerical Instability! Determinant value is too small: " + determinant[kkk] +
                        "Try running with fewer annotations and then with fewer Gaussians.");
            }
        }
        final double mult[] = new double[numAnnotations];
        double sumWeight = 0.0;
        for( int iii = 0; iii < data.length; iii++ ) {
            sumWeight += data[iii].weight;
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
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

                pVarInClusterLog10[kkk] = pClusterLog10[kkk] + (( -0.5 * sum )/Math.log(10.0)) - denomLog10[kkk];
                final double pVar = Math.pow(10.0, pVarInClusterLog10[kkk]);
                likelihood += pVar * data[iii].weight;

                if( pVarInClusterLog10[kkk] > 0.0 || Double.isNaN(pVarInClusterLog10[kkk]) || Double.isInfinite(pVarInClusterLog10[kkk]) ) {
                    logger.warn("det = " + sigma[kkk].det());
                    logger.warn("denom = " + denomLog10[kkk]);
                    logger.warn("sumExp = " + sum);
                    logger.warn("pVar = " + pVar);
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                            logger.warn(sigmaVals[kkk][ppp][jjj]);
                        }
                    }

                    logger.warn("About to throw exception due to numerical instability. Try running with fewer annotations and then with fewer Gaussians. " +
                            "It is best to only use the annotations which appear to be Gaussianly distributed for this Gaussian mixture model.");
                    throw new StingException("Numerical Instability! Found NaN after performing log10: " + pVarInClusterLog10[kkk] + ", cluster = " + kkk + ", variant index = " + iii);
                }
            }

            pVarInClusterReals = MathUtils.normalizeFromLog10( pVarInClusterLog10 );
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                pVarInCluster[kkk][iii] = pVarInClusterReals[kkk] * data[iii].weight;
                if( Double.isNaN(pVarInCluster[kkk][iii]) ) {
                    logger.warn("About to throw exception due to numerical instability. Try running with fewer annotations and then with fewer Gaussians. " +
                            "It is best to only use the annotations which appear to be Gaussianly distributed for this Gaussian mixture model.");
                    throw new StingException("Numerical Instability! Found NaN after rescaling log10 values: " + pVarInCluster[kkk][iii] + ", cluster = " + kkk + ", variant index = " + iii);
                }
            }
        }

        logger.info("Explained likelihood = " + String.format("%.5f",likelihood / sumWeight));
        return likelihood / sumWeight;
    }


    private void evaluateGaussiansForSingleVariant( final double[] annotations, final double[] pVarInCluster ) {

        final int numAnnotations = annotations.length;
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

            final double denomLog10 = Math.log10(Math.pow(2.0 * Math.PI, ((double)numAnnotations) / 2.0)) + Math.log10(Math.pow(determinant[kkk], 0.5));
            pVarInCluster[kkk] = Math.pow(10.0, pClusterLog10[kkk] + (( -0.5 * sum ) / Math.log(10.0)) - denomLog10);
        }
    }


    private void maximizeGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;
        final double sigmaVals[][][] = new double[numGaussians][numAnnotations][numAnnotations];
        final double meanVals[][] = new double[numGaussians][numAnnotations];

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                meanVals[kkk][jjj] = 0.0;
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] = 0.0;
                }
            }
        }
        double sumPK = 0.0;
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            double sumProb = 0.0;
            for( int iii = 0; iii < numVariants; iii++ ) {
                final double prob = pVarInCluster[kkk][iii];
                if( Double.isNaN(prob) ) {
                    logger.warn("About to throw exception due to numerical instability. Try running with fewer annotations and then with fewer Gaussians. " +
                            "It is best to only use the annotations which appear to be Gaussianly distributed for this Gaussian mixture model.");
                    throw new StingException("Numerical Instability! Found NaN in M-step: " + pVarInCluster[kkk][iii] + ", cluster = " + kkk + ", variant index = " + iii);
                }
                sumProb += prob;
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    meanVals[kkk][jjj] +=  prob * data[iii].annotations[jjj];
                }
            }

            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                meanVals[kkk][jjj] /=  sumProb;
            }

            for( int iii = 0; iii < numVariants; iii++ ) {
                final double prob = pVarInCluster[kkk][iii];
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                        sigmaVals[kkk][jjj][ppp] +=  prob * (data[iii].annotations[jjj]-meanVals[kkk][jjj]) * (data[iii].annotations[ppp]-meanVals[kkk][ppp]);
                    }
                }
            }

            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                    if( sigmaVals[kkk][jjj][ppp] < MIN_SIGMA && sigmaVals[kkk][jjj][ppp] > -MIN_SIGMA ) { // Very small numbers are a very big problem
                        logger.warn("The sigma values look exceptionally small.... Probably about to crash due to numeric instability.");
                    }
                    sigmaVals[kkk][ppp][jjj] = sigmaVals[kkk][jjj][ppp]; // sigma must be a symmetric matrix
                }
            }

            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] /= sumProb;
                }
            }

            if( FORCE_INDEPENDENT_ANNOTATIONS ) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        if(jjj!=ppp) {
                            sigmaVals[kkk][jjj][ppp] = 0.0;
                        }
                    }
                }

            }
            
            final Matrix tmpMatrix = new Matrix(sigmaVals[kkk]);
            if( tmpMatrix.det() > MIN_DETERMINANT ) {
                sigma[kkk] = new Matrix(sigmaVals[kkk]);
                determinant[kkk] = sigma[kkk].det();

                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    mu[kkk][jjj] = meanVals[kkk][jjj];
                }
            } else {
                logger.warn("Tried to create a covariance matrix with exceptionally small determinant.");
            }

            pClusterLog10[kkk] = sumProb;
            sumPK += sumProb;
        }

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
           pClusterLog10[kkk] = Math.log10( pClusterLog10[kkk] / sumPK );
        }

        pClusterLog10 = MathUtils.normalizeFromLog10( pClusterLog10, true );
    }
}