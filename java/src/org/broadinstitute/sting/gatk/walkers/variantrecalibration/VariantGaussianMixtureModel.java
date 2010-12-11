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
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import Jama.*; 

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public final class VariantGaussianMixtureModel extends VariantOptimizationModel {

    protected final static Logger logger = Logger.getLogger(VariantGaussianMixtureModel.class);

    public final VariantDataManager dataManager;
    private final int maxGaussians;
    private final int maxIterations;
    private final static long RANDOM_SEED = 91801305;
    private final Random rand = new Random( RANDOM_SEED );
    private final double MIN_PROB_CONVERGENCE = 1E-5;

    private final double SHRINKAGE;
    private final double DIRICHLET_PARAMETER;
    private final boolean FORCE_INDEPENDENT_ANNOTATIONS;

    private final double[][] mu; // The means for each cluster
    private final Matrix[] sigma; // The covariance matrix for each cluster
    private final double[][][] sigmaInverse;
    private double[] pClusterLog10;
    private final double[] determinant;
    private final double[] sqrtDeterminantLog10;
    private final double stdThreshold;
    private double singletonFPRate = -1; // Estimated FP rate for singleton calls.  Used to estimate FP rate as a function of AC

    private double[] empiricalMu;
    private Matrix empiricalSigma;

    private final double[] hyperParameter_a;
    private final double[] hyperParameter_b;
    private final double[] hyperParameter_lambda;

    private final double CONSTANT_GAUSSIAN_DENOM_LOG10;

    private static final Pattern COMMENT_PATTERN = Pattern.compile("^##.*");
    private static final Pattern ANNOTATION_PATTERN = Pattern.compile("^@!ANNOTATION.*");
    private static final Pattern CLUSTER_PATTERN = Pattern.compile("^@!CLUSTER.*");

    public VariantGaussianMixtureModel() {
        dataManager = null;
        maxGaussians = 0;
        maxIterations = 0;
        SHRINKAGE = 0;
        DIRICHLET_PARAMETER = 0;
        FORCE_INDEPENDENT_ANNOTATIONS = false;
        mu = null;
        sigma = null;
        sigmaInverse = null;
        determinant = null;
        sqrtDeterminantLog10 = null;
        stdThreshold = 0;
        hyperParameter_a = null;
        hyperParameter_b = null;
        hyperParameter_lambda = null;
        CONSTANT_GAUSSIAN_DENOM_LOG10 = 0.0;
    }

    public VariantGaussianMixtureModel( final VariantDataManager _dataManager, final int _maxGaussians, final int _maxIterations,
                                        final boolean _forceIndependent, final double _stdThreshold, final double _shrinkage, final double _dirichlet) {
        dataManager = _dataManager;
        maxGaussians = _maxGaussians;
        maxIterations = _maxIterations;

        mu = new double[maxGaussians][];
        sigma = new Matrix[maxGaussians];
        determinant = new double[maxGaussians];
        sqrtDeterminantLog10 = null;
        pClusterLog10 = new double[maxGaussians];
        stdThreshold = _stdThreshold;
        FORCE_INDEPENDENT_ANNOTATIONS = _forceIndependent;
        hyperParameter_a = new double[maxGaussians];
        hyperParameter_b = new double[maxGaussians];
        hyperParameter_lambda = new double[maxGaussians];
        sigmaInverse = null; // This field isn't used during GenerateVariantClusters pass

        CONSTANT_GAUSSIAN_DENOM_LOG10 =  Math.log10(Math.pow(2.0 * Math.PI, ((double)dataManager.numAnnotations) / 2.0));
        SHRINKAGE = _shrinkage;
        DIRICHLET_PARAMETER = _dirichlet;
    }

    public VariantGaussianMixtureModel( final double _targetTITV, final File clusterFile, final double backOffGaussianFactor ) {
        super( _targetTITV );
        final ExpandingArrayList<String> annotationLines = new ExpandingArrayList<String>();
        final ExpandingArrayList<String> clusterLines = new ExpandingArrayList<String>();

        try {
            for ( final String line : new XReadLines( clusterFile ) ) {
                if( ANNOTATION_PATTERN.matcher(line).matches() ) {
                    annotationLines.add(line);
                } else if( CLUSTER_PATTERN.matcher(line).matches() ) {
                    clusterLines.add(line);
                } else if( !COMMENT_PATTERN.matcher(line).matches() ) {
                    throw new UserException.MalformedFile(clusterFile, "Could not parse line: " + line);
                }
            }
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(clusterFile, e);
        }

        dataManager = new VariantDataManager( annotationLines );
        // Several of the clustering parameters aren't used the second time around in VariantRecalibrator.java
        SHRINKAGE = 0;
        DIRICHLET_PARAMETER = 0;
        maxIterations = 0;
        stdThreshold = 0.0;
        FORCE_INDEPENDENT_ANNOTATIONS = false;
        hyperParameter_a = null;
        hyperParameter_b = null;
        hyperParameter_lambda = null;
        determinant = null;

        // BUGBUG: move this parsing out of the constructor
        CONSTANT_GAUSSIAN_DENOM_LOG10 =  Math.log10(Math.pow(2.0 * Math.PI, ((double)dataManager.numAnnotations) / 2.0));
        maxGaussians = clusterLines.size();
        mu = new double[maxGaussians][dataManager.numAnnotations];
        final double sigmaVals[][][] = new double[maxGaussians][dataManager.numAnnotations][dataManager.numAnnotations];
        sigma = new Matrix[maxGaussians];
        sigmaInverse = new double[maxGaussians][dataManager.numAnnotations][dataManager.numAnnotations];
        pClusterLog10 = new double[maxGaussians];
        sqrtDeterminantLog10 = new double[maxGaussians];

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
            sigmaInverse[kkk] = sigma[kkk].inverse().getArray(); // Precompute all the inverses and determinants for use later
            sqrtDeterminantLog10[kkk] = Math.log10(Math.pow(sigma[kkk].det(), 0.5)); // Precompute for use later
            kkk++;
        }

        logger.info("Found " + maxGaussians + " clusters using " + dataManager.numAnnotations + " annotations: " + dataManager.annotationKeys);
    }
    
    public final void run( final PrintStream clusterFile ) {

        // Only cluster with a good set of knowns. Filter based on being too many std's away from the mean annotation value
        // Filtering based on known status and qual threshold happens in GenerateVariantClusters
        long numAboveSTD = 0L;
        for( int iii = 0; iii < dataManager.data.length; iii++ ) {
            final VariantDatum datum = dataManager.data[iii];
            for( final double val : datum.annotations ) {
                if( Math.abs(val) > stdThreshold ) {
                    datum.weight = 0.0;
                    numAboveSTD++;
                    break;
                }
            }
        }

        logger.info( numAboveSTD + " variants were rejected from the training set for having annotation values more than X standard deviations away from the mean. (--stdThreshold = " + stdThreshold + ")" );

        generateEmpricalStats( dataManager.data );

        logger.info("Initializing using k-means...");
        initializeUsingKMeans( dataManager.data );
        logger.info("... done!");
        createClusters( dataManager.data, 0, maxGaussians, clusterFile );
    }

    private void generateEmpricalStats( final VariantDatum[] data ) {
        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        empiricalMu = new double[numAnnotations];
        final double[][] sigmaVals = new double[numAnnotations][numAnnotations];


        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            empiricalMu[jjj] = 0.0;
            for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                sigmaVals[jjj][ppp] = 0.0;
            }
        }

        for(int iii = 0; iii < numVariants; iii++) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                empiricalMu[jjj] += data[iii].annotations[jjj] / ((double) numVariants);
            }
        }

        for(int iii = 0; iii < numVariants; iii++) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    if( jjj == ppp ) {
                        sigmaVals[jjj][ppp] = 1.0; //0.01 * numAnnotations;
                    } else {
                        sigmaVals[jjj][ppp] = 0.0;
                    }
                    //sigmaVals[jjj][ppp] += (data[iii].annotations[jjj]-empiricalMu[jjj]) * (data[iii].annotations[ppp]-empiricalMu[ppp]);
                }
            }
        }

        empiricalSigma = new Matrix(sigmaVals);
        //empiricalSigma.timesEquals(1.0 / (Math.pow(maxGaussians, 2.0 / ((double) numAnnotations))));
    }

    public void setSingletonFPRate( final double rate ) {
        this.singletonFPRate = rate;
    }

    public final double getAlleleCountPrior( final int alleleCount, final double maxProbability ) {
        return Math.min( maxProbability, 1.0 - Math.pow(singletonFPRate, alleleCount) );
    }

    private void initializeUsingKMeans( final VariantDatum[] data ) {
        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        for( int kkk = 0; kkk < maxGaussians; kkk++ ) {
            mu[kkk] = data[rand.nextInt(numVariants)].annotations;
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mu[kkk][jjj] = -4.0 + 8.0 * rand.nextDouble();
            }
        }

        for( int ttt = 0; ttt < 60; ttt++ ) {
            performKMeansIteration( data );
        }
    }

    private void performKMeansIteration( final VariantDatum[] data ) {
        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final int[] assignment = new int[numVariants];
        for( int iii = 0; iii < numVariants; iii++ ) {
            final VariantDatum datum = data[iii];
            if(datum.weight > 0.0) {
                double minDistance = Double.MAX_VALUE;
                int minCluster = -1;
                for( int kkk = 0; kkk < maxGaussians; kkk++ ) {
                    double dist = 0.0;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        dist += (datum.annotations[jjj] - mu[kkk][jjj]) * (datum.annotations[jjj] - mu[kkk][jjj]);
                    }
                    if(dist < minDistance) {
                        minDistance = dist;
                        minCluster = kkk;
                    }
                }
                assignment[iii] = minCluster;
            }
        }

        for( int kkk = 0; kkk < maxGaussians; kkk++ ) {
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mu[kkk][jjj] = 0.0;
            }
            int numAssigned = 0;

            for( int iii = 0; iii < numVariants; iii++ ) {
                if(assignment[iii] == kkk) {
                    numAssigned++;
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        mu[kkk][jjj] += data[iii].annotations[jjj];
                    }
                }
            }
            if(numAssigned != 0) {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    mu[kkk][jjj] /= (double) numAssigned;
                }
            } else {
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                   mu[kkk][jjj] = -4.0 + 8.0 * rand.nextDouble();
                }
            }
        }
    }

    public final void createClusters( final VariantDatum[] data, final int startCluster, final int stopCluster, final PrintStream clusterFile ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;

        final double[][] pVarInCluster = new double[maxGaussians][numVariants]; // Probability that the variant is in that cluster = simply evaluate the multivariate Gaussian

        // loop control variables:
        // iii - loop over data points
        // jjj - loop over annotations (features)
        // ppp - loop over annotations again (full rank covariance matrix)
        // kkk - loop over clusters
        // ttt - loop over EM iterations

        // Set up the initial random Gaussians
        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            hyperParameter_a[kkk] = numAnnotations;
            hyperParameter_b[kkk] = SHRINKAGE;
            hyperParameter_lambda[kkk] = DIRICHLET_PARAMETER;
            pClusterLog10[kkk] = Math.log10(1.0 / ((double) (stopCluster - startCluster)));
            final double[][] randSigma = new double[numAnnotations][numAnnotations];
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = jjj; ppp < numAnnotations; ppp++ ) {
                    randSigma[ppp][jjj] = 0.55 + 1.25 * rand.nextDouble();
                    if(rand.nextBoolean()) {
                        randSigma[ppp][jjj] *= -1.0;
                    }
                    if(jjj != ppp) { randSigma[jjj][ppp] = 0.0; } // Sigma is a symmetric, positive-definite matrix created by taking a lower diagonal matrix and multiplying it by its transpose
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

            logger.info("Finished iteration " + ttt );
            ttt++;
            if( Math.abs(currentLikelihood - previousLikelihood) < MIN_PROB_CONVERGENCE ) {
                logger.info("Convergence!");
                break;
            }
            previousLikelihood = currentLikelihood;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Output the final cluster parameters
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printClusterParameters( clusterFile );
    }

    private void printClusterParameters( final PrintStream clusterFile ) {
        dataManager.printClusterFileHeader( clusterFile );

        final int numAnnotations = mu[0].length;
        for( int kkk = 0; kkk < maxGaussians; kkk++ ) {
            if( Math.pow(10.0, pClusterLog10[kkk]) > 1E-4 ) { // BUGBUG: make this a command line argument
                final double sigmaVals[][] = sigma[kkk].getArray();
                clusterFile.print("@!CLUSTER");
                clusterFile.print(String.format(",%.8f", Math.pow(10.0, pClusterLog10[kkk])));
                for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    clusterFile.print(String.format(",%.8f", mu[kkk][jjj]));
                }
                for(int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for(int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        clusterFile.print(String.format(",%.8f", (sigmaVals[jjj][ppp] / hyperParameter_a[kkk]) ));
                    }
                }
                clusterFile.println();
            }
        }
    }

    public double decodeAnnotation(GenomeLocParser genomeLocParser, final String annotationKey, final VariantContext vc, final boolean jitter ) {
        double value;
        if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // HRun values must be jittered a bit to work in this GMM
            value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            value += -0.25 + 0.5 * rand.nextDouble();
        } else if( annotationKey.equals("QUAL") ) {
            value = vc.getPhredScaledQual();
        } else {
            try {
                value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            } catch( Exception e ) {
                throw new UserException.MalformedFile(vc.getSource(), "No double value detected for annotation = " + annotationKey + " in variant at " + VariantContextUtils.getLocation(genomeLocParser,vc) + ", reported annotation value = " + vc.getAttribute( annotationKey ), e );
            }
        }
        return value;
    }

    public final double evaluateVariant( GenomeLocParser genomeLocParser, final VariantContext vc ) {
        final double[] log10pVarInCluster = new double[maxGaussians];
        final double[] annotations = new double[dataManager.numAnnotations];

        for( int jjj = 0; jjj < dataManager.numAnnotations; jjj++ ) {
            final double value = decodeAnnotation( genomeLocParser, dataManager.annotationKeys.get(jjj), vc, false );
            annotations[jjj] = (value - dataManager.meanVector[jjj]) / dataManager.varianceVector[jjj];
        }

        return evaluateGaussiansForSingleVariant( annotations, log10pVarInCluster );
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Code to determine FDR tranches for VariantDatum[]
    //
    // ---------------------------------------------------------------------------------------------------------
    public static abstract class SelectionMetric {
        String name = null;

        public SelectionMetric(String name) {
            this.name = name;
        }

        public String getName() { return name; }

        public abstract double getThreshold(double tranche);
        public abstract double getTarget();
        public abstract void calculateRunningMetric(List<VariantDatum> data);
        public abstract double getRunningMetric(int i);
        public abstract int datumValue(VariantDatum d);
    }

    public static class NovelTiTvMetric extends SelectionMetric {
        double[] runningTiTv;
        double targetTiTv = 0;

        public NovelTiTvMetric(double target) {
            super("NovelTiTv");
            targetTiTv = target; // compute the desired TiTv
        }

        public double getThreshold(double tranche) {
            return fdrToTiTv(tranche, targetTiTv);
        }

        public double getTarget() { return targetTiTv; }

        public void calculateRunningMetric(List<VariantDatum> data) {
            int ti = 0, tv = 0;
            runningTiTv = new double[data.size()];

            for ( int i = data.size() - 1; i >= 0; i-- ) {
                VariantDatum datum = data.get(i);
                if ( ! datum.isKnown ) {
                    if ( datum.isTransition ) { ti++; } else { tv++; }
                    runningTiTv[i] = ti / Math.max(1.0 * tv, 1.0);
                }
            }
        }

        public double getRunningMetric(int i) {
            return runningTiTv[i];
        }

        public int datumValue(VariantDatum d) {
            return d.isTransition ? 1 : 0;
        }
    }

    public static class TruthSensitivityMetric extends SelectionMetric {
        double[] runningSensitivity;
        double targetTiTv = 0;
        int nTrueSites = 0;

        public TruthSensitivityMetric(int nTrueSites) {
            super("TruthSensitivity");
            this.nTrueSites = nTrueSites;
        }

        public double getThreshold(double tranche) {
            return tranche/100; // tranche of 1 => 99% sensivity target
        }

        public double getTarget() { return 1.0; }

        public void calculateRunningMetric(List<VariantDatum> data) {
            int nCalledAtTruth = 0;
            runningSensitivity = new double[data.size()];

            for ( int i = data.size() - 1; i >= 0; i-- ) {
                VariantDatum datum = data.get(i);
                nCalledAtTruth += datum.atTruthSite ? 1 : 0;
                runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
            }
        }

        public double getRunningMetric(int i) {
            return runningSensitivity[i];
        }

        public int datumValue(VariantDatum d) {
            return d.atTruthSite ? 1 : 0;
        }
    }

    public final static List<Tranche> findTranches( final VariantDatum[] data, final double[] tranches, SelectionMetric metric ) {
        return findTranches( data, tranches, metric, null );
    }

    public final static List<Tranche> findTranches( final VariantDatum[] data, final double[] trancheThresholds, SelectionMetric metric, File debugFile ) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.length, data.length));

        List<VariantDatum> tranchesData = sortVariantsbyLod(data);
        metric.calculateRunningMetric(tranchesData);

        if ( debugFile != null) { writeTranchesDebuggingInfo(debugFile, tranchesData, metric); }

        List<Tranche> tranches = new ArrayList<Tranche>();
        for ( double trancheThreshold : trancheThresholds ) {
            Tranche t = findTranche(tranchesData, metric, trancheThreshold);

            if ( t == null ) {
                if ( tranches.size() == 0 )
                    throw new UserException(String.format("Couldn't find any tranche containing variants with a %s > %.2f", metric.getName(), metric.getThreshold(trancheThreshold)));
                break;
            }

            tranches.add(t);
        }

        return tranches;
    }

    private static void writeTranchesDebuggingInfo(File f, List<VariantDatum> tranchesData, SelectionMetric metric ) {
        try {
            PrintStream out = new PrintStream(f);
            out.println("Qual metricValue runningValue");
            for ( int i = 0; i < tranchesData.size(); i++ ) {
                VariantDatum  d = tranchesData.get(i);
                int score = metric.datumValue(d);
                double runningValue = metric.getRunningMetric(i);
                out.printf("%.4f %d %.4f%n", d.lod, score, runningValue);
            }
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(f, e);
        }
    }

    private final static List<VariantDatum> sortVariantsbyLod(final VariantDatum[] data) {
        List<VariantDatum> sorted = new ArrayList<VariantDatum>(Arrays.asList(data));
        Collections.sort(sorted);
        return sorted;
    }

    public final static Tranche findTranche( final List<VariantDatum> data, final SelectionMetric metric, final double trancheThreshold ) {
        logger.info(String.format("  Tranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        double metricThreshold = metric.getThreshold(trancheThreshold);
        int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( metric.getRunningMetric(i) >= metricThreshold ) {
                // we've found the largest group of variants with Ti/Tv >= our target titv
                Tranche t = trancheOfVariants(data, i, trancheThreshold, metric.getTarget());
                logger.info(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.info(String.format("  Tranche is %s", t));
                return t;
            }
        }

        // we get here when there's no subset of variants with Ti/Tv >= threshold, in which case we should return null
        //logger.info(String.format("  **Couldn't find tranche for %.3f FDR = %.3f TiTv threshold; last running TiTv was %.2f ",
        //        desiredFDR, titvThreshold, runningTiTv[runningTiTv.length-1]));


//        for ( int i = 0; i < runningTiTv.length; i++ ) {
//            if ( runningTiTv[i] >= titvThreshold ) {
//                // we've found the largest group of variants with Ti/Tv >= our target titv
//                Tranche t = trancheOfVariants(data, i, desiredFDR, targetTiTv);
//                logger.info(String.format("  Found tranche for %.3f FDR = %.3f TiTv threshold starting with variant %d; running TiTv is %.2f ",
//                        desiredFDR, titvThreshold, i, metric.getRunningMetric(i)));
//                logger.info(String.format("  Tranche is %s", t));
//                return t;
//            }
//        }
//
//        // we get here when there's no subset of variants with Ti/Tv >= threshold, in which case we should return null
//        logger.info(String.format("  **Couldn't find tranche for %.3f FDR = %.3f TiTv threshold; last running TiTv was %.2f ",
//                desiredFDR, titvThreshold, runningTiTv[runningTiTv.length-1]));
//
        return null;
    }

    public final static Tranche trancheOfVariants( final List<VariantDatum> data, int minI, double fdr, double target ) {
        int numKnown = 0, numNovel = 0, knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;

        double minLod = data.get(minI).lod;
        VariantDatum last = null;
        for ( VariantDatum datum : data ) {
            if ( datum.lod >= minLod ) {
                //if( ! datum.isKnown ) System.out.println(datum.pos);
                if ( datum.isKnown ) {
                    numKnown++;
                    if ( datum.isTransition ) { knownTi++; } else { knownTv++; }
                } else {
                    numNovel++;
                    if ( datum.isTransition ) { novelTi++; } else { novelTv++; }

                }
            }
            last = datum;
        }

        double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
        double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

        int accessibleTruthSites = countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        int nCallsAtTruth = countCallsAtTruth(data, minLod);


        return new Tranche(fdr, target, minLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, nCallsAtTruth);
    }

    public final static double fdrToTiTv(double desiredFDR, double targetTiTv) {
            return (1.0 - desiredFDR / 100.0) * (targetTiTv - 0.5) + 0.5;
    }

    public static int countCallsAtTruth(final VariantDatum[] data, double minLOD ) {
        int n = 0;
        for ( VariantDatum d : data) { n += d.atTruthSite && d.lod >= minLOD ? 1 : 0; }
        return n;
    }

    public static int countCallsAtTruth(final List<VariantDatum> data, double minLOD ) {
        int n = 0;
        for ( VariantDatum d : data) { n += d.atTruthSite && d.lod >= minLOD ? 1 : 0; }
        return n;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Code to evaluate vectors given an already trained model
    //
    // ---------------------------------------------------------------------------------------------------------
    private double evaluateGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numAnnotations = data[0].annotations.length;
        double likelihood = 0.0;
        final double sigmaVals[][][] = new double[maxGaussians][][];
        final double denomLog10[] = new double[maxGaussians];
        final double pVarInClusterLog10[] = new double[maxGaussians];
        double pVarInClusterReals[];

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            sigmaVals[kkk] = sigma[kkk].inverse().getArray();
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    sigmaVals[kkk][jjj][ppp] *= hyperParameter_a[kkk];
                }
            }
            double sum = 0.0;
            for(int jjj = 1; jjj < numAnnotations; jjj++) {
                sum += diGamma((hyperParameter_a[kkk] + 1.0 - jjj) / 2.0);
            }
            sum -= Math.log(determinant[kkk]);
            sum += Math.log(2.0) * numAnnotations;
            final double gamma = 0.5 * sum;
            sum = 0.0;
            for(int ccc = 0; ccc < maxGaussians; ccc++) {
                sum += hyperParameter_lambda[ccc];
            }
            final double pi = diGamma(hyperParameter_lambda[kkk]) - diGamma(sum);
            final double beta = (-1.0 * numAnnotations) / (2.0 * hyperParameter_b[kkk]);
            denomLog10[kkk] = (pi / Math.log(10.0)) + (gamma / Math.log(10.0)) + (beta / Math.log(10.0));
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

                pVarInClusterLog10[kkk] = (( -0.5 * sum )/Math.log(10.0)) + denomLog10[kkk];
                final double pVar = Math.pow(10.0, pVarInClusterLog10[kkk]);
                likelihood += pVar * data[iii].weight;

                if( Double.isNaN(pVarInClusterLog10[kkk]) || Double.isInfinite(pVarInClusterLog10[kkk]) ) {
                    logger.warn("det = " + sigma[kkk].det());
                    logger.warn("denom = " + denomLog10[kkk]);
                    logger.warn("sumExp = " + sum);
                    logger.warn("mixtureLog10 = " + pClusterLog10[kkk]);
                    logger.warn("pVar = " + pVar);
                    for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                        for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                            logger.warn(sigmaVals[kkk][ppp][jjj]);
                        }
                    }

                    logger.warn("About to throw exception due to numerical instability. Try running with fewer annotations and then with fewer Gaussians. " +
                            "It is best to only use the annotations which appear to be Gaussianly distributed for this Gaussian mixture model.");
                    throw new ReviewedStingException("Numerical Instability! Found NaN after performing log10: " + pVarInClusterLog10[kkk] + ", cluster = " + kkk + ", variant index = " + iii);
                }
            }

            pVarInClusterReals = MathUtils.normalizeFromLog10( pVarInClusterLog10 );
            for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
                pVarInCluster[kkk][iii] = pVarInClusterReals[kkk] * data[iii].weight;
                if( Double.isNaN(pVarInCluster[kkk][iii]) ) {
                    logger.warn("About to throw exception due to numerical instability. Try running with fewer annotations and then with fewer Gaussians. " +
                            "It is best to only use the annotations which appear to be Gaussianly distributed for this Gaussian mixture model.");
                    throw new ReviewedStingException("Numerical Instability! Found NaN after rescaling log10 values: " + pVarInCluster[kkk][iii] + ", cluster = " + kkk + ", variant index = " + iii);
                }
            }
        }

        logger.info("explained likelihood = " + String.format("%.5f",likelihood / sumWeight));
        return likelihood / sumWeight;
    }

    private double evaluateGaussiansForSingleVariant( final double[] annotations, final double[] log10pVarInCluster ) {

        final int numAnnotations = annotations.length;

        for( int kkk = 0; kkk < maxGaussians; kkk++ ) {
            final double sigmaVals[][] = sigmaInverse[kkk];
            double sum = 0.0;
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                double value = 0.0;
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    final double myMu = mu[kkk][ppp];
                    final double myAnn = annotations[ppp];
                    final double mySigma = sigmaVals[ppp][jjj];
                    value += (myAnn - myMu) * mySigma;
                }
                final double jNorm = annotations[jjj] - mu[kkk][jjj];
                final double prod = value * jNorm;
                sum += prod;
            }

            final double denomLog10 = CONSTANT_GAUSSIAN_DENOM_LOG10 + sqrtDeterminantLog10[kkk];
            final double evalGaussianPDFLog10 = (( -0.5 * sum ) / Math.log(10.0)) - denomLog10;
            //final double pVar1 = Math.pow(10.0, pClusterLog10[kkk] + evalGaussianPDFLog10);
            final double pVar1 = pClusterLog10[kkk] + evalGaussianPDFLog10;
            log10pVarInCluster[kkk] = pVar1;
        }

        return MathUtils.log10sumLog10(log10pVarInCluster);
    }


    private void maximizeGaussians( final VariantDatum[] data, final double[][] pVarInCluster, final int startCluster, final int stopCluster ) {

        final int numVariants = data.length;
        final int numAnnotations = data[0].annotations.length;
        final double sigmaVals[][][] = new double[maxGaussians][numAnnotations][numAnnotations];
        final double wishartVals[][] = new double[numAnnotations][numAnnotations];
        final double meanVals[][] = new double[maxGaussians][numAnnotations];

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
                    throw new ReviewedStingException("Numerical Instability! Found NaN in M-step: " + pVarInCluster[kkk][iii] + ", cluster = " + kkk + ", variant index = " + iii);
                }
                sumProb += prob;
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    meanVals[kkk][jjj] += prob * data[iii].annotations[jjj];
                }
            }

            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                meanVals[kkk][jjj] = (meanVals[kkk][jjj] + SHRINKAGE * empiricalMu[jjj]) / (sumProb + SHRINKAGE);
            }

            final double shrinkageFactor = (SHRINKAGE * sumProb) / (SHRINKAGE + sumProb);
            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                    wishartVals[jjj][ppp] = shrinkageFactor * (meanVals[kkk][jjj] - empiricalMu[jjj]) * (meanVals[kkk][ppp] - empiricalMu[ppp]);
                }
            }

            for( int iii = 0; iii < numVariants; iii++ ) {
                final double prob = pVarInCluster[kkk][iii];
                for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                    for( int ppp = 0; ppp < numAnnotations; ppp++ ) {
                        sigmaVals[kkk][jjj][ppp] += prob * (data[iii].annotations[jjj]-meanVals[kkk][jjj]) * (data[iii].annotations[ppp]-meanVals[kkk][ppp]);
                    }
                }
            }

            final Matrix tmpMatrix = empiricalSigma.plus(new Matrix(wishartVals).plus(new Matrix(sigmaVals[kkk])));

            sigma[kkk] = (Matrix)tmpMatrix.clone();
            determinant[kkk] = sigma[kkk].det();

            for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
                mu[kkk][jjj] = meanVals[kkk][jjj];
            }

            pClusterLog10[kkk] = sumProb;
            sumPK += sumProb;

            hyperParameter_a[kkk] = sumProb + numAnnotations;
            hyperParameter_b[kkk] = sumProb + SHRINKAGE;
            hyperParameter_lambda[kkk] = sumProb + DIRICHLET_PARAMETER;
        }

        for( int kkk = startCluster; kkk < stopCluster; kkk++ ) {
            pClusterLog10[kkk] = Math.log10( pClusterLog10[kkk] / sumPK );
        }

        pClusterLog10 = MathUtils.normalizeFromLog10( pClusterLog10, true );
    }

    // from http://en.wikipedia.org/wiki/Digamma_function
    // According to J.M. Bernardo AS 103 algorithm the digamma function for x, a real number, can be approximated by:
    private static double diGamma(final double x) {
        return Math.log(x) - ( 1.0 / (2.0 * x) )
                           - ( 1.0 / (12.0 * Math.pow(x, 2.0)) )
                           + ( 1.0 / (120.0 * Math.pow(x, 4.0)) )
                           - ( 1.0 / (252.0 * Math.pow(x, 6.0)) ); 
    }
}
