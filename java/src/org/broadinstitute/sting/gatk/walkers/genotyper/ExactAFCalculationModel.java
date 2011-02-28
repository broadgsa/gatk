/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;
import java.io.PrintStream;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {
    private final static boolean DEBUG = false;
    private final static boolean PRINT_LIKELIHOODS = false;

    public enum ExactCalculation {
        N2_GOLD_STANDARD,
        LINEAR_EXPERIMENTAL
    }

    private final static boolean COMPARE_TO_GS = true;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    private boolean SIMPLE_GREEDY_GENOTYPER = false;
    private static final double[] log10Cache;
    private static final double[] jacobianLogTable;
    private static final int JACOBIAN_LOG_TABLE_SIZE = 101;
    private static final double JACOBIAN_LOG_TABLE_STEP = 0.1;
    private static final double MAX_JACOBIAN_TOLERANCE = 10.0;
    private static final int MAXN = 10000;

    static {
        log10Cache = new double[2*MAXN];
        jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k=1; k < 2*MAXN; k++)
            log10Cache[k] = Math.log10(k);

        for (int k=0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
            jacobianLogTable[k] = Math.log10(1.0+Math.pow(10.0,-((double)k)
                    * JACOBIAN_LOG_TABLE_STEP));

        }

    }

    final private ExactCalculation calcToUse;
    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        calcToUse = UAC.EXACT_CALCULATION_TYPE;
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, Genotype> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors) {
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        double[] gsPosteriors;
        if ( COMPARE_TO_GS ) // due to annoying special values in incoming array, we have to clone up here
            gsPosteriors = log10AlleleFrequencyPosteriors.clone();

        int lastK = -1;
        switch ( calcToUse ) {
            case N2_GOLD_STANDARD:
                lastK = gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);
                break;
            case LINEAR_EXPERIMENTAL:
                lastK = linearExact(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);
                break;
        }

        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        if ( COMPARE_TO_GS ) {
            gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, gsPosteriors);

            double log10thisPVar = Math.log10(MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors)[0]);
            double log10gsPVar = Math.log10(MathUtils.normalizeFromLog10(gsPosteriors)[0]);
            boolean eq = (log10thisPVar == Double.NEGATIVE_INFINITY && log10gsPVar == Double.NEGATIVE_INFINITY) || MathUtils.compareDoubles(log10thisPVar, log10gsPVar, 1e-4) == 0;

            if ( ! eq || PRINT_LIKELIHOODS ) {
                System.out.printf("----------------------------------------%n");
                for (int k=0; k < log10AlleleFrequencyPosteriors.length; k++) {
                    System.out.printf("  %d\t%.2f\t%.2f\t%b%n", k,
                            log10AlleleFrequencyPosteriors[k], gsPosteriors[k],
                            log10AlleleFrequencyPosteriors[k] == gsPosteriors[k]);
                }
                System.out.printf("MAD_AC\t%d\t%d\t%.2f\t%.2f\t%.6f%n",
                        ref.getLocus().getStart(), lastK, log10thisPVar, log10gsPVar, log10thisPVar - log10gsPVar);
            }
        }

    }

    private static final double[][] getGLs(Map<String, Genotype> GLs) {
        double[][] genotypeLikelihoods = new double[GLs.size()+1][];

        int j = 0;
        for ( Genotype sample : GLs.values() ) {
            j++;

            if ( sample.hasLikelihoods() ) {
                //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
                genotypeLikelihoods[j] = sample.getLikelihoods().getAsVector();
            }
        }

        return genotypeLikelihoods;
    }

    private static class ExactACCache {
        double[] kMinus2, kMinus1, kMinus0;

        private static double[] create(int n, double defaultValue) {
            double[] v = new double[n];
            Arrays.fill(v, defaultValue);
            return v;
        }

        public ExactACCache(int n, double defaultValue) {
            kMinus2 = create(n, defaultValue);
            kMinus1 = create(n, defaultValue);
            kMinus0 = create(n, defaultValue);
        }

        public void rotate() {
            double[] tmp = kMinus2;
            kMinus2 = kMinus1;
            kMinus1 = kMinus0;
            kMinus0 = tmp;
        }

        public double[] getkMinus2() {
            return kMinus2;
        }

        public double[] getkMinus1() {
            return kMinus1;
        }

        public double[] getkMinus0() {
            return kMinus0;
        }
    }

    public int linearExact(Map<String, Genotype> GLs,
                           double[] log10AlleleFrequencyPriors,
                           double[] log10AlleleFrequencyPosteriors) {
        int numSamples = GLs.size();
        int numChr = 2*numSamples;
        double[][] genotypeLikelihoods = getGLs(GLs); // todo -- remove me, not sure this is helping
        double[] logNumerator = new double[3];

        // set posteriors to negative infinity by default:
        //Arrays.fill(log10AlleleFrequencyPosteriors, Double.NEGATIVE_INFINITY);

        ExactACCache logY = new ExactACCache(numSamples+1, Double.NEGATIVE_INFINITY);
        logY.getkMinus0()[0] = 0.0; // the zero case

        double maxLog10L = Double.NEGATIVE_INFINITY;
        boolean done = false;
        int lastK = -1;

        // todo -- we may be able to start second loop some way down the calculation, since GdAs loop only
        // todo -- considers part of the matrix as well
        for (int k=0; k <= numChr && ! done; k++ ) {
            double[] kMinus0 = logY.getkMinus0();
            double[] kMinus1 = logY.getkMinus1();
            double[] kMinus2 = logY.getkMinus2();

            if ( k == 0 ) {
                // special case for k = 0
                for ( int j=1; j <= numSamples; j++ ) {
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods[j][GenotypeType.AA.ordinal()];
                }
            } else { // k > 0
                for ( int j=1; j <= numSamples; j++ ) {
                    double[] gl = genotypeLikelihoods[j];
                    double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];

                    if (k < 2*j-1)
                        logNumerator[0] = log10Cache[2*j-k] + log10Cache[2*j-k-1] + kMinus0[j-1] +
                                gl[GenotypeType.AA.ordinal()];
                    else
                        logNumerator[0] = Double.NEGATIVE_INFINITY;

                    if (k < 2*j)
                        logNumerator[1] = log10Cache[2*k] + log10Cache[2*j-k]+ kMinus1[j-1] +
                                gl[GenotypeType.AB.ordinal()];
                    else
                        logNumerator[1] = Double.NEGATIVE_INFINITY;

                    if (k > 1)
                        logNumerator[2] = log10Cache[k] + log10Cache[k-1] + kMinus2[j-1] +
                                gl[GenotypeType.BB.ordinal()];
                    else
                        logNumerator[2] = Double.NEGATIVE_INFINITY;

                    kMinus0[j] = softMax(logNumerator) - logDenominator;
                }
            }

            // update the posteriors vector
            double log10LofK = kMinus0[numSamples];
            log10AlleleFrequencyPosteriors[k] = log10LofK + log10AlleleFrequencyPriors[k];

            // can we abort early?
            lastK = k;
            maxLog10L = Math.max(maxLog10L, log10LofK);
            if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
                if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
                done = true;
            }

            logY.rotate();
        }

        return lastK;
    }

    public int gdaN2GoldStandard(Map<String, Genotype> GLs,
                                 double[] log10AlleleFrequencyPriors,
                                 double[] log10AlleleFrequencyPosteriors) {
        int numSamples = GLs.size();
        int numChr = 2*numSamples;

        double[][] logYMatrix = new double[1+numSamples][1+numChr];

        for (int i=0; i <=numSamples; i++)
            for (int j=0; j <=numChr; j++)
                logYMatrix[i][j] = Double.NEGATIVE_INFINITY;

        //YMatrix[0][0] = 1.0;
        logYMatrix[0][0] = 0.0;
        int j=0;

        for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
            j++;

            if ( !sample.getValue().hasLikelihoods() )
                continue;

            //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
            double[] genotypeLikelihoods = sample.getValue().getLikelihoods().getAsVector();
            //double logDenominator = Math.log10(2.0*j*(2.0*j-1));
            double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];

            // special treatment for k=0: iteration reduces to:
            //YMatrix[j][0] = YMatrix[j-1][0]*genotypeLikelihoods[GenotypeType.AA.ordinal()];
            logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[GenotypeType.AA.ordinal()];

            for (int k=1; k <= 2*j; k++ ) {

                //double num = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[GenotypeType.AA.ordinal()];
                double logNumerator[];
                logNumerator = new double[3];
                if (k < 2*j-1)
                    logNumerator[0] = log10Cache[2*j-k] + log10Cache[2*j-k-1] + logYMatrix[j-1][k] +
                            genotypeLikelihoods[GenotypeType.AA.ordinal()];
                else
                    logNumerator[0] = Double.NEGATIVE_INFINITY;


                if (k < 2*j)
                    logNumerator[1] = log10Cache[2*k] + log10Cache[2*j-k]+ logYMatrix[j-1][k-1] +
                            genotypeLikelihoods[GenotypeType.AB.ordinal()];
                else
                    logNumerator[1] = Double.NEGATIVE_INFINITY;

                if (k > 1)
                    logNumerator[2] = log10Cache[k] + log10Cache[k-1] + logYMatrix[j-1][k-2] +
                            genotypeLikelihoods[GenotypeType.BB.ordinal()];
                else
                    logNumerator[2] = Double.NEGATIVE_INFINITY;

                double logNum = softMax(logNumerator);

                //YMatrix[j][k] = num/den;
                logYMatrix[j][k] = logNum - logDenominator;
            }

        }

        for (int k=0; k <= numChr; k++)
            log10AlleleFrequencyPosteriors[k] = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];

        return numChr;
    }

    private final static void printLikelihoods(int numChr, double[][] logYMatrix, double[] log10AlleleFrequencyPriors) {
        int j = logYMatrix.length - 1;
        System.out.printf("-----------------------------------%n");
        for (int k=0; k <= numChr; k++) {
            double posterior = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];
            System.out.printf("  %4d\t%8.2f\t%8.2f\t%8.2f%n", k, logYMatrix[j][k], log10AlleleFrequencyPriors[k], posterior);
        }
    }

    double softMax(double[] vec) {
        // compute naively log10(10^x[0] + 10^x[1]+...)
        //        return Math.log10(MathUtils.sumLog10(vec));

        // better approximation: do Jacobian logarithm function on data pairs
        double a = softMaxPair(vec[0],vec[1]);
        return softMaxPair(a,vec[2]);
    }

    static public double softMaxPair(double x, double y) {
        if (Double.isInfinite(x))
            return y;

        if (Double.isInfinite(y))
            return x;

        if (y >= x + MAX_JACOBIAN_TOLERANCE)
            return y;
        if (x >= y + MAX_JACOBIAN_TOLERANCE)
            return x;

        // OK, so |y-x| < tol: we use the following identity then:
        // we need to compute log10(10^x + 10^y)
        // By Jacobian logarithm identity, this is equal to
        // max(x,y) + log10(1+10^-abs(x-y))
        // we compute the second term as a table lookup
        // with integer quantization
        double diff = Math.abs(x-y);
        double t1 =x;
        if (y > x)
            t1 = y;
        // t has max(x,y)
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        int ind = (int)Math.round(diff/JACOBIAN_LOG_TABLE_STEP);
        double t2 = jacobianLogTable[ind];

        // gdebug+
        //double z =Math.log10(1+Math.pow(10.0,-diff));
        //System.out.format("x: %f, y:%f, app: %f, true: %f ind:%d\n",x,y,t2,z,ind);
        //gdebug-
        return t1+t2;
        // return Math.log10(Math.pow(10.0,x) + Math.pow(10.0,y));
    }



    /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    public Map<String, Genotype> assignGenotypes(VariantContext vc,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int AFofMaxLikelihood) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        Allele refAllele = vc.getReference();
        Allele altAllele = vc.getAlternateAllele(0);

        Map<String, Genotype> GLs = vc.getGenotypes();
        double[][] pathMetricArray = new double[GLs.size()+1][AFofMaxLikelihood+1];
        int[][] tracebackArray = new int[GLs.size()+1][AFofMaxLikelihood+1];

        ArrayList<String> sampleIndices = new ArrayList<String>();
        int sampleIdx = 0;

        // todo - optimize initialization
        for (int k=0; k <= AFofMaxLikelihood; k++)
            for (int j=0; j <= GLs.size(); j++)
                pathMetricArray[j][k] = -1e30;

        pathMetricArray[0][0] = 0.0;

        if (SIMPLE_GREEDY_GENOTYPER) {
            sampleIndices.addAll(GLs.keySet());
            sampleIdx = GLs.size();
        }
        else {

            for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
                if ( !sample.getValue().hasLikelihoods() )
                    continue;

                double[] likelihoods = sample.getValue().getLikelihoods().getAsVector();
                sampleIndices.add(sample.getKey());

                for (int k=0; k <= AFofMaxLikelihood; k++) {

                    double bestMetric = pathMetricArray[sampleIdx][k] + likelihoods[0];
                    int bestIndex = k;

                    if (k>0) {
                        double m2 =  pathMetricArray[sampleIdx][k-1] + likelihoods[1];
                        if (m2 > bestMetric) {
                            bestMetric = m2;
                            bestIndex  = k-1;
                        }
                    }

                    if (k>1) {
                        double m2 =  pathMetricArray[sampleIdx][k-2] + likelihoods[2];
                        if (m2 > bestMetric) {
                            bestMetric = m2;
                            bestIndex  = k-2;
                        }
                    }

                    pathMetricArray[sampleIdx+1][k] = bestMetric;
                    tracebackArray[sampleIdx+1][k] = bestIndex;
                }
                sampleIdx++;
            }
        }

        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();

        int startIdx = AFofMaxLikelihood;
        for (int k = sampleIdx; k > 0; k--) {
            int bestGTguess;
            String sample = sampleIndices.get(k-1);
            Genotype g = GLs.get(sample);
            if ( !g.hasLikelihoods() )
                continue;

            if (SIMPLE_GREEDY_GENOTYPER)
                bestGTguess = Utils.findIndexOfMaxEntry(g.getLikelihoods().getAsVector());
            else {
                int newIdx = tracebackArray[k][startIdx];
                bestGTguess = startIdx - newIdx;
                startIdx = newIdx;
            }

            ArrayList<Allele> myAlleles = new ArrayList<Allele>();

            double qual;
            double[] likelihoods = g.getLikelihoods().getAsVector();

            if (bestGTguess == 0) {
                myAlleles.add(refAllele);
                myAlleles.add(refAllele);
                qual = likelihoods[0] - Math.max(likelihoods[1], likelihoods[2]);
            } else if(bestGTguess == 1) {
                myAlleles.add(refAllele);
                myAlleles.add(altAllele);
                qual = likelihoods[1] - Math.max(likelihoods[0], likelihoods[2]);

            }  else {
                myAlleles.add(altAllele);
                myAlleles.add(altAllele);
                qual = likelihoods[2] - Math.max(likelihoods[1], likelihoods[0]);
            }


            if (qual < 0) {
                // QUAL can be negative if the chosen genotype is not the most likely one individually.
                // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen on
                double[] normalized = MathUtils.normalizeFromLog10(likelihoods);
                double chosenGenotype = normalized[bestGTguess];
                qual = -1.0 * Math.log10(1.0 - chosenGenotype);
            }

            calls.put(sample, new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));

        }

        return calls;
    }

}

// working linearized version
//public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {
//    private final static boolean PRINT_LIKELIHOODS = false;
//
//    public enum ExactCalculation {
//        N2_GOLD_STANDARD,
//        LINEAR_EXPERIMENTAL
//    }
//
//    private final static boolean COMPARE_TO_GS = false;
//    private final static boolean PRINT_MAD_AC_POSTERIORS = false;
//    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
//
//    private boolean SIMPLE_GREEDY_GENOTYPER = false;
//    private static final double[] log10Cache;
//    private static final double[] jacobianLogTable;
//    private static final int JACOBIAN_LOG_TABLE_SIZE = 101;
//    private static final double JACOBIAN_LOG_TABLE_STEP = 0.1;
//    private static final double MAX_JACOBIAN_TOLERANCE = 10.0;
//    private static final int MAXN = 10000;
//
//    static {
//        log10Cache = new double[2*MAXN];
//	    jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];
//
//        log10Cache[0] = Double.NEGATIVE_INFINITY;
//        for (int k=1; k < 2*MAXN; k++)
//            log10Cache[k] = Math.log10(k);
//
//	for (int k=0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
//	    jacobianLogTable[k] = Math.log10(1.0+Math.pow(10.0,-((double)k)
//						       * JACOBIAN_LOG_TABLE_STEP));
//
//	}
//
//    }
//
//    final private ExactCalculation calcToUse;
//    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
//        super(UAC, N, logger, verboseWriter);
//        calcToUse = UAC.EXACT_CALCULATION_TYPE;
//    }
//
//    public void getLog10PNonRef(RefMetaDataTracker tracker,
//                                ReferenceContext ref,
//                                Map<String, Genotype> GLs,
//                                double[] log10AlleleFrequencyPriors,
//                                double[] log10AlleleFrequencyPosteriors) {
//        switch ( calcToUse ) {
//            case N2_GOLD_STANDARD:
//                gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);
//                break;
//            case LINEAR_EXPERIMENTAL:
//                madByAC(ref, GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);
//                break;
//        }
//    }
//
//    private static final double[][] getGLs(Map<String, Genotype> GLs) {
//        double[][] genotypeLikelihoods = new double[GLs.size()+1][];
//
//        int j = 0;
//        for ( Genotype sample : GLs.values() ) {
//            j++;
//
//            if ( sample.hasLikelihoods() ) {
//                //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
//                genotypeLikelihoods[j] = sample.getLikelihoods().getAsVector();
//            }
//        }
//
//        return genotypeLikelihoods;
//    }
//
//    private static class ExactACCache {
//        double[] kMinus2, kMinus1, kMinus0;
//
//        private static double[] create(int n, double defaultValue) {
//            double[] v = new double[n];
//            Arrays.fill(v, defaultValue);
//            return v;
//        }
//
//        public ExactACCache(int nSamples, double defaultValue) {
//            kMinus2 = create(nSamples, defaultValue);
//            kMinus1 = create(nSamples, defaultValue);
//            kMinus0 = create(nSamples, defaultValue);
//        }
//
//        public void rotate() {
//            double[] tmp = kMinus2;
//            kMinus2 = kMinus1;
//            kMinus1 = kMinus0;
//            kMinus0 = tmp;
//        }
//
//        public double[] getkMinus2() {
//            return kMinus2;
//        }
//
//        public double[] getkMinus1() {
//            return kMinus1;
//        }
//
//        public double[] getkMinus0() {
//            return kMinus0;
//        }
//    }
//
//    public void madByAC(ReferenceContext ref,
//                        Map<String, Genotype> GLs,
//                        double[] log10AlleleFrequencyPriors,
//                        double[] log10AlleleFrequencyPosteriors) {
//        // todo -- remove me after testing
//        double[] gsPosteriors = log10AlleleFrequencyPosteriors;
//        if ( COMPARE_TO_GS ) {
//            gsPosteriors = log10AlleleFrequencyPosteriors.clone();
//            gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, gsPosteriors);
//        }
//
//        int numSamples = GLs.size();
//        int numChr = 2*numSamples;
//        double[][] genotypeLikelihoods = getGLs(GLs); // todo -- remove me, not sure this is helping
//
//        // set posteriors to negative infinity by default:
//        Arrays.fill(log10AlleleFrequencyPosteriors, Double.NEGATIVE_INFINITY);
//
//        // todo -- replace this matrix with 3 vectors (k, k-1, k-2) and cycle through them
//        // todo -- this is *CRITICAL* to reduce the algorithm to a true ~linear algorithm
//        double[][] logYMatrix = new double[1+numSamples][1+numChr];
//        for (int i=0; i <= numSamples; i++) // initialize
//            Arrays.fill(logYMatrix[i], Double.NEGATIVE_INFINITY);
//        logYMatrix[0][0] = 0.0; // the zero case
//
//        double maxLog10L = Double.NEGATIVE_INFINITY;
//        boolean done = false;
//        int lastK = -1;
//
//        // todo -- we may be able to start second loop some way down the calculation, since GdAs loop only
//        // todo -- considers part of the matrix as well
//        for (int k=0; k <= numChr && ! done; k++ ) {
//            if ( k == 0 ) {
//                // special case for k = 0
//                for ( int j=1; j <= numSamples; j++ ) {
//                    logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[j][GenotypeType.AA.ordinal()];
//                }
//            } else { // k > 0
//                for ( int j=1; j <= numSamples; j++ ) {
//                    double[] gl = genotypeLikelihoods[j];
//                    double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];
//
//                    double[] logNumerator = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
//                    if (k < 2*j-1)
//                        logNumerator[0] = log10Cache[2*j-k] + log10Cache[2*j-k-1] + logYMatrix[j-1][k] +
//                                gl[GenotypeType.AA.ordinal()];
//
//                    if (k < 2*j)
//                        logNumerator[1] = log10Cache[2*k] + log10Cache[2*j-k]+ logYMatrix[j-1][k-1] +
//                                gl[GenotypeType.AB.ordinal()];
//
//                    if (k > 1)
//                        logNumerator[2] = log10Cache[k] + log10Cache[k-1] + logYMatrix[j-1][k-2] +
//                                gl[GenotypeType.BB.ordinal()];
//
//                    logYMatrix[j][k] = softMax(logNumerator) - logDenominator;
//                }
//            }
//
//            // update the posteriors vector
//            double log10LofK = logYMatrix[numSamples][k];
//            log10AlleleFrequencyPosteriors[k] = log10LofK + log10AlleleFrequencyPriors[k];
//
//            // can we abort early?
//            lastK = k;
//            maxLog10L = Math.max(maxLog10L, log10LofK);
//            if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
//                if ( PRINT_MAD_AC_POSTERIORS )
//                    System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
//                done = true;
//            }
//        }
//
//        if ( PRINT_MAD_AC_POSTERIORS ) {
//            System.out.printf("----------------------------------------%n");
//            for (int k=0; k <= numChr; k++) {
//                System.out.printf("  %d\t%.2f\t%.2f\t%b%n", k,
//                        log10AlleleFrequencyPosteriors[k], gsPosteriors[k],
//                        log10AlleleFrequencyPosteriors[k] == gsPosteriors[k]);
//            }
//            double log10thisPVar = Math.log10(MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors)[0]);
//            double log10gsPVar = Math.log10(MathUtils.normalizeFromLog10(gsPosteriors)[0]);
//            System.out.printf("MAD_AC\t%d\t%d\t%.2f\t%.2f\t%.6f%n",
//                    ref.getLocus().getStart(), lastK, log10thisPVar, log10gsPVar, log10thisPVar - log10gsPVar);
//        }
//
//        if ( PRINT_LIKELIHOODS ) printLikelihoods(numChr, logYMatrix, log10AlleleFrequencyPriors);
//    }
//
//
//    public void gdaN2GoldStandard(Map<String, Genotype> GLs,
//                                  double[] log10AlleleFrequencyPriors,
//                                  double[] log10AlleleFrequencyPosteriors) {
//        int numSamples = GLs.size();
//        int numChr = 2*numSamples;
//
//        double[][] logYMatrix = new double[1+numSamples][1+numChr];
//
//        for (int i=0; i <=numSamples; i++)
//            for (int j=0; j <=numChr; j++)
//                logYMatrix[i][j] = Double.NEGATIVE_INFINITY;
//
//        //YMatrix[0][0] = 1.0;
//        logYMatrix[0][0] = 0.0;
//        int j=0;
//
//        for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
//            j++;
//
//            if ( !sample.getValue().hasLikelihoods() )
//                continue;
//
//            //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
//            double[] genotypeLikelihoods = sample.getValue().getLikelihoods().getAsVector();
//            //double logDenominator = Math.log10(2.0*j*(2.0*j-1));
//            double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];
//
//            // special treatment for k=0: iteration reduces to:
//            //YMatrix[j][0] = YMatrix[j-1][0]*genotypeLikelihoods[GenotypeType.AA.ordinal()];
//            logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[GenotypeType.AA.ordinal()];
//
//            for (int k=1; k <= 2*j; k++ ) {
//
//                //double num = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[GenotypeType.AA.ordinal()];
//                double logNumerator[];
//                logNumerator = new double[3];
//                if (k < 2*j-1)
//                    logNumerator[0] = log10Cache[2*j-k] + log10Cache[2*j-k-1] + logYMatrix[j-1][k] +
//                        genotypeLikelihoods[GenotypeType.AA.ordinal()];
//                else
//                    logNumerator[0] = Double.NEGATIVE_INFINITY;
//
//
//                if (k < 2*j)
//                    logNumerator[1] = log10Cache[2*k] + log10Cache[2*j-k]+ logYMatrix[j-1][k-1] +
//                        genotypeLikelihoods[GenotypeType.AB.ordinal()];
//                else
//                    logNumerator[1] = Double.NEGATIVE_INFINITY;
//
//                if (k > 1)
//                    logNumerator[2] = log10Cache[k] + log10Cache[k-1] + logYMatrix[j-1][k-2] +
//                        genotypeLikelihoods[GenotypeType.BB.ordinal()];
//                else
//                    logNumerator[2] = Double.NEGATIVE_INFINITY;
//
//                double logNum = softMax(logNumerator);
//
//                //YMatrix[j][k] = num/den;
//                logYMatrix[j][k] = logNum - logDenominator;
//            }
//
//        }
//
//
//        for (int k=0; k <= numChr; k++)
//            log10AlleleFrequencyPosteriors[k] = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];
//
//        if ( PRINT_LIKELIHOODS ) printLikelihoods(numChr, logYMatrix, log10AlleleFrequencyPriors);
//    }
//
//    private final static void printLikelihoods(int numChr, double[][] logYMatrix, double[] log10AlleleFrequencyPriors) {
//        int j = logYMatrix.length - 1;
//        System.out.printf("-----------------------------------%n");
//        for (int k=0; k <= numChr; k++) {
//            double posterior = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];
//            System.out.printf("  %4d\t%8.2f\t%8.2f\t%8.2f%n", k, logYMatrix[j][k], log10AlleleFrequencyPriors[k], posterior);
//        }
//    }
//
//    double softMax(double[] vec) {
//        // compute naively log10(10^x[0] + 10^x[1]+...)
//        //        return Math.log10(MathUtils.sumLog10(vec));
//
//        // better approximation: do Jacobian logarithm function on data pairs
//        double a = softMaxPair(vec[0],vec[1]);
//        return softMaxPair(a,vec[2]);
//    }
//
//    static public double softMaxPair(double x, double y) {
//        if (Double.isInfinite(x))
//            return y;
//
//        if (Double.isInfinite(y))
//            return x;
//
//        if (y >= x + MAX_JACOBIAN_TOLERANCE)
//            return y;
//        if (x >= y + MAX_JACOBIAN_TOLERANCE)
//            return x;
//
//        // OK, so |y-x| < tol: we use the following identity then:
//        // we need to compute log10(10^x + 10^y)
//        // By Jacobian logarithm identity, this is equal to
//        // max(x,y) + log10(1+10^-abs(x-y))
//        // we compute the second term as a table lookup
//        // with integer quantization
//        double diff = Math.abs(x-y);
//        double t1 =x;
//        if (y > x)
//            t1 = y;
//        // t has max(x,y)
//        // we have pre-stored correction for 0,0.1,0.2,... 10.0
//        int ind = (int)Math.round(diff/JACOBIAN_LOG_TABLE_STEP);
//        double t2 = jacobianLogTable[ind];
//
//        // gdebug+
//        //double z =Math.log10(1+Math.pow(10.0,-diff));
//        //System.out.format("x: %f, y:%f, app: %f, true: %f ind:%d\n",x,y,t2,z,ind);
//        //gdebug-
//        return t1+t2;
//        // return Math.log10(Math.pow(10.0,x) + Math.pow(10.0,y));
//    }
//
//
//
//    /**
//     * Can be overridden by concrete subclasses
//     * @param vc                   variant context with genotype likelihoods
//     * @param log10AlleleFrequencyPosteriors    allele frequency results
//     * @param AFofMaxLikelihood    allele frequency of max likelihood
//     *
//     * @return calls
//     */
//    public Map<String, Genotype> assignGenotypes(VariantContext vc,
//                                                 double[] log10AlleleFrequencyPosteriors,
//                                                 int AFofMaxLikelihood) {
//        if ( !vc.isVariant() )
//            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());
//
//        Allele refAllele = vc.getReference();
//        Allele altAllele = vc.getAlternateAllele(0);
//
//        Map<String, Genotype> GLs = vc.getGenotypes();
//        double[][] pathMetricArray = new double[GLs.size()+1][AFofMaxLikelihood+1];
//        int[][] tracebackArray = new int[GLs.size()+1][AFofMaxLikelihood+1];
//
//        ArrayList<String> sampleIndices = new ArrayList<String>();
//        int sampleIdx = 0;
//
//        // todo - optimize initialization
//        for (int k=0; k <= AFofMaxLikelihood; k++)
//            for (int j=0; j <= GLs.size(); j++)
//                pathMetricArray[j][k] = -1e30;
//
//        pathMetricArray[0][0] = 0.0;
//
//        if (SIMPLE_GREEDY_GENOTYPER) {
//            sampleIndices.addAll(GLs.keySet());
//            sampleIdx = GLs.size();
//        }
//        else {
//
//            for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
//                if ( !sample.getValue().hasLikelihoods() )
//                    continue;
//
//                double[] likelihoods = sample.getValue().getLikelihoods().getAsVector();
//                sampleIndices.add(sample.getKey());
//
//                for (int k=0; k <= AFofMaxLikelihood; k++) {
//
//                    double bestMetric = pathMetricArray[sampleIdx][k] + likelihoods[0];
//                    int bestIndex = k;
//
//                    if (k>0) {
//                        double m2 =  pathMetricArray[sampleIdx][k-1] + likelihoods[1];
//                        if (m2 > bestMetric) {
//                            bestMetric = m2;
//                            bestIndex  = k-1;
//                        }
//                    }
//
//                    if (k>1) {
//                        double m2 =  pathMetricArray[sampleIdx][k-2] + likelihoods[2];
//                        if (m2 > bestMetric) {
//                            bestMetric = m2;
//                            bestIndex  = k-2;
//                        }
//                    }
//
//                    pathMetricArray[sampleIdx+1][k] = bestMetric;
//                    tracebackArray[sampleIdx+1][k] = bestIndex;
//                }
//                sampleIdx++;
//            }
//        }
//
//        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();
//
//        int startIdx = AFofMaxLikelihood;
//        for (int k = sampleIdx; k > 0; k--) {
//            int bestGTguess;
//            String sample = sampleIndices.get(k-1);
//            Genotype g = GLs.get(sample);
//            if ( !g.hasLikelihoods() )
//                continue;
//
//            if (SIMPLE_GREEDY_GENOTYPER)
//                bestGTguess = Utils.findIndexOfMaxEntry(g.getLikelihoods().getAsVector());
//            else {
//                int newIdx = tracebackArray[k][startIdx];
//                bestGTguess = startIdx - newIdx;
//                startIdx = newIdx;
//            }
//
//            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
//
//            double qual;
//            double[] likelihoods = g.getLikelihoods().getAsVector();
//
//            if (bestGTguess == 0) {
//                myAlleles.add(refAllele);
//                myAlleles.add(refAllele);
//                qual = likelihoods[0] - Math.max(likelihoods[1], likelihoods[2]);
//            } else if(bestGTguess == 1) {
//                myAlleles.add(refAllele);
//                myAlleles.add(altAllele);
//                qual = likelihoods[1] - Math.max(likelihoods[0], likelihoods[2]);
//
//            }  else {
//                myAlleles.add(altAllele);
//                myAlleles.add(altAllele);
//                qual = likelihoods[2] - Math.max(likelihoods[1], likelihoods[0]);
//            }
//
//
//            if (qual < 0) {
//                // QUAL can be negative if the chosen genotype is not the most likely one individually.
//                // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen on
//                double[] normalized = MathUtils.normalizeFromLog10(likelihoods);
//                double chosenGenotype = normalized[bestGTguess];
//                qual = -1.0 * Math.log10(1.0 - chosenGenotype);
//            }
//
//            calls.put(sample, new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));
//
//        }
//
//        return calls;
//    }
//
//}