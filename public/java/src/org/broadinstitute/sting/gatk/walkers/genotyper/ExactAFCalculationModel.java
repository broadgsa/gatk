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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {
    //
    // code for testing purposes
    //
    private final static boolean DEBUG = false;
    private final static boolean PRINT_LIKELIHOODS = false;
    private final static int N_CYCLES = 1;
    private SimpleTimer timerExpt = new SimpleTimer("linearExactBanded");
    private SimpleTimer timerGS = new SimpleTimer("linearExactGS");
    private final static boolean COMPARE_TO_GS = false;

    public enum ExactCalculation {
        N2_GOLD_STANDARD,
        LINEAR_EXPERIMENTAL
    }

    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    private boolean SIMPLE_GREEDY_GENOTYPER = false;

    private final static double SUM_GL_THRESH_NOCALL = -0.001; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    final private ExactCalculation calcToUse;
    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        calcToUse = UAC.EXACT_CALCULATION_TYPE;
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, Genotype> GLs, Set<Allele>alleles,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors) {
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        double[] gsPosteriors;
        if ( COMPARE_TO_GS ) // due to annoying special values in incoming array, we have to clone up here
            gsPosteriors = log10AlleleFrequencyPosteriors.clone();

        int idxAA = GenotypeType.AA.ordinal();
        int idxAB = GenotypeType.AB.ordinal();
        int idxBB = GenotypeType.BB.ordinal();

        // todo -- remove me after testing
        if ( N_CYCLES > 1 ) {
            for ( int i = 0; i < N_CYCLES; i++) {
                timerGS.restart();
                linearExact(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.clone(), idxAA, idxAB, idxBB);
                timerGS.stop();

                timerExpt.restart();
                linearExactBanded(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.clone());
                timerExpt.stop();
            }

            System.out.printf("good = %.2f, expt = %.2f, delta = %.2f%n",
                    timerGS.getElapsedTime(), timerExpt.getElapsedTime(), timerExpt.getElapsedTime()-timerGS.getElapsedTime());
        }

        int lastK = -1;

        int numAlleles = alleles.size();

        int idxDiag = numAlleles;
        int incr = numAlleles - 1;

        double[][] posteriorCache = new double[numAlleles-1][];
        double[] bestAFguess = new double[numAlleles-1];

        for (int k=1; k < numAlleles; k++) {
            // multi-allelic approximation, part 1: Ideally
            // for each alt allele compute marginal (suboptimal) posteriors -
            // compute indices for AA,AB,BB for current allele - genotype likelihoods are a linear vector that can be thought of
            // as a row-wise upper triangular matrix of likelihoods.
            // So, for example, with 2 alt alleles, likelihoods have AA,AB,AC,BB,BC,CC.
            // 3 alt alleles: AA,AB,AC,AD BB BC BD CC CD DD

            idxAA = 0;
            idxAB = k;
            // yy is always element on the diagonal.
            // 2 alleles: BBelement 2
            // 3 alleles: BB element  3. CC element 5
            // 4 alleles:
            idxBB = idxDiag;
            idxDiag += incr--;

            // todo - possible cleanup
            switch ( calcToUse ) {
                case N2_GOLD_STANDARD:
                    lastK = gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors, idxAA, idxAB, idxBB);
                    break;
                case LINEAR_EXPERIMENTAL:
                    lastK = linearExact(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors, idxAA, idxAB, idxBB);
                    break;
            }
            if (numAlleles > 2) {
                posteriorCache[k-1] = log10AlleleFrequencyPosteriors.clone();
                bestAFguess[k-1] = (double)MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors);
            }
        }

        if (numAlleles > 2) {
            // multiallelic approximation, part 2:
            // report posteriors for allele that has highest estimated AC
            int mostLikelyAlleleIdx = MathUtils.maxElementIndex(bestAFguess);
            for (int k=0; k < log10AlleleFrequencyPosteriors.length-1; k++)
                log10AlleleFrequencyPosteriors[k] = (posteriorCache[mostLikelyAlleleIdx][k]);

        }
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        // todo -- REMOVE ME AFTER TESTING
        if ( COMPARE_TO_GS ) {
            gdaN2GoldStandard(GLs, log10AlleleFrequencyPriors, gsPosteriors, idxAA, idxAB, idxBB);

            double log10thisPVar = Math.log10(MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors)[0]);
            double log10gsPVar = Math.log10(MathUtils.normalizeFromLog10(gsPosteriors)[0]);
            boolean eq = (log10thisPVar == Double.NEGATIVE_INFINITY && log10gsPVar == Double.NEGATIVE_INFINITY) || MathUtils.compareDoubles(log10thisPVar, log10gsPVar, 1e-4) == 0;

            if ( ! eq || PRINT_LIKELIHOODS ) {
                System.out.printf("----------------------------------------%n");
                for (int k=0; k < log10AlleleFrequencyPosteriors.length; k++) {
                    double x = log10AlleleFrequencyPosteriors[k];
                    System.out.printf("  %d\t%.2f\t%.2f\t%b%n", k,
                            x < -1e10 ? Double.NEGATIVE_INFINITY : x, gsPosteriors[k],
                            log10AlleleFrequencyPosteriors[k] == gsPosteriors[k]);
                }
                System.out.printf("MAD_AC\t%d\t%d\t%.2f\t%.2f\t%.6f%n",
                        ref.getLocus().getStart(), lastK, log10thisPVar, log10gsPVar, log10thisPVar - log10gsPVar);
            }
        }

    }

    private static final ArrayList<double[]> getGLs(Map<String, Genotype> GLs) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>();

        //int j = 0;
        genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.values() ) {
            if ( sample.hasLikelihoods() ) {
                //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
                double[] gls = sample.getLikelihoods().getAsVector();

                if (MathUtils.sum(gls) < SUM_GL_THRESH_NOCALL)
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }


    // -------------------------------------------------------------------------------------
    //
    // Linearized, ~O(N), implementation.
    //
    // -------------------------------------------------------------------------------------

    /**
     * A simple data structure that holds the current, prev, and prev->prev likelihoods vectors
     * for the exact model calculation
     */
    private final static class ExactACCache {
        double[] kMinus2, kMinus1, kMinus0;

        private final static double[] create(int n) {
            return new double[n];
        }

        public ExactACCache(int n) {
            kMinus2 = create(n);
            kMinus1 = create(n);
            kMinus0 = create(n);
        }

        final public void rotate() {
            double[] tmp = kMinus2;
            kMinus2 = kMinus1;
            kMinus1 = kMinus0;
            kMinus0 = tmp;
        }

        final public double[] getkMinus2() {
            return kMinus2;
        }

        final public double[] getkMinus1() {
            return kMinus1;
        }

        final public double[] getkMinus0() {
            return kMinus0;
        }
    }

    // now with banding
    public int linearExactBanded(Map<String, Genotype> GLs,
                                 double[] log10AlleleFrequencyPriors,
                                 double[] log10AlleleFrequencyPosteriors) {
        throw new NotImplementedException();
//        final int numSamples = GLs.size();
//        final int numChr = 2*numSamples;
//        final double[][] genotypeLikelihoods = getGLs(GLs);
//
//        final ExactACCache logY = new ExactACCache(numSamples+1);
//        logY.getkMinus0()[0] = 0.0; // the zero case
//
//        double maxLog10L = Double.NEGATIVE_INFINITY;
//        boolean done = false;
//        int lastK = -1;
//        final int BAND_SIZE = 10;
//
//        for (int k=0; k <= numChr && ! done; k++ ) {
//            final double[] kMinus0 = logY.getkMinus0();
//            int jStart = Math.max(k - BAND_SIZE, 1);
//            int jStop = Math.min(k + BAND_SIZE, numSamples);
//
//            if ( k == 0 ) { // special case for k = 0
//                for ( int j=1; j <= numSamples; j++ ) {
//                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods[j][GenotypeType.AA.ordinal()];
//                }
//            } else { // k > 0
//                final double[] kMinus1 = logY.getkMinus1();
//                final double[] kMinus2 = logY.getkMinus2();
//                Arrays.fill(kMinus0,0);
//
//                for ( int j = jStart; j <= jStop; j++ ) {
//                    final double[] gl = genotypeLikelihoods[j];
//                    final double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];
//
//                    double aa = Double.NEGATIVE_INFINITY;
//                    double ab = Double.NEGATIVE_INFINITY;
//                    if (k < 2*j-1)
//                        aa = log10Cache[2*j-k] + log10Cache[2*j-k-1] + kMinus0[j-1] + gl[GenotypeType.AA.ordinal()];
//
//                    if (k < 2*j)
//                        ab = log10Cache[2*k] + log10Cache[2*j-k]+ kMinus1[j-1] + gl[GenotypeType.AB.ordinal()];
//
//                    double log10Max;
//                    if (k > 1) {
//                        final double bb = log10Cache[k] + log10Cache[k-1] + kMinus2[j-1] + gl[GenotypeType.BB.ordinal()];
//                        log10Max = approximateLog10SumLog10(aa, ab, bb);
//                    } else {
//                        // we know we aren't considering the BB case, so we can use an optimized log10 function
//                        log10Max = approximateLog10SumLog10(aa, ab);
//                    }
//
//                    // finally, update the L(j,k) value
//                    kMinus0[j] = log10Max - logDenominator;
//
//                    String offset = Utils.dupString(' ',k);
//                    System.out.printf("%s%3d %3d %.2f%n", offset, k, j, kMinus0[j]);
//                }
//            }
//
//            // update the posteriors vector
//            final double log10LofK = kMinus0[jStop];
//            log10AlleleFrequencyPosteriors[k] = log10LofK + log10AlleleFrequencyPriors[k];
//
//            // can we abort early?
//            lastK = k;
//            maxLog10L = Math.max(maxLog10L, log10LofK);
//            if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
//                if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
//                done = true;
//            }
//
//            logY.rotate();
//        }
//
//        return lastK;
    }

    public int linearExact(Map<String, Genotype> GLs,
                           double[] log10AlleleFrequencyPriors,
                           double[] log10AlleleFrequencyPosteriors, int idxAA, int idxAB, int idxBB) {
        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        final ExactACCache logY = new ExactACCache(numSamples+1);
        logY.getkMinus0()[0] = 0.0; // the zero case

        double maxLog10L = Double.NEGATIVE_INFINITY;
        boolean done = false;
        int lastK = -1;

        for (int k=0; k <= numChr && ! done; k++ ) {
            final double[] kMinus0 = logY.getkMinus0();

            if ( k == 0 ) { // special case for k = 0
                for ( int j=1; j <= numSamples; j++ ) {
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods.get(j)[idxAA];
                }
            } else { // k > 0
                final double[] kMinus1 = logY.getkMinus1();
                final double[] kMinus2 = logY.getkMinus2();

                for ( int j=1; j <= numSamples; j++ ) {
                    final double[] gl = genotypeLikelihoods.get(j);
                    final double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];

                    double aa = Double.NEGATIVE_INFINITY;
                    double ab = Double.NEGATIVE_INFINITY;
                    if (k < 2*j-1)
                        aa = MathUtils.log10Cache[2*j-k] + MathUtils.log10Cache[2*j-k-1] + kMinus0[j-1] + gl[idxAA];

                    if (k < 2*j)
                        ab = MathUtils.log10Cache[2*k] + MathUtils.log10Cache[2*j-k]+ kMinus1[j-1] + gl[idxAB];

                    double log10Max;
                    if (k > 1) {
                        final double bb = MathUtils.log10Cache[k] + MathUtils.log10Cache[k-1] + kMinus2[j-1] + gl[idxBB];
                        log10Max = approximateLog10SumLog10(aa, ab, bb);
                    } else {
                        // we know we aren't considering the BB case, so we can use an optimized log10 function
                        log10Max = approximateLog10SumLog10(aa, ab);
                    }

                    // finally, update the L(j,k) value
                    kMinus0[j] = log10Max - logDenominator;
                }
            }

            // update the posteriors vector
            final double log10LofK = kMinus0[numSamples];
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

    final static double approximateLog10SumLog10(double a, double b, double c) {
        //return softMax(new double[]{a, b, c});
        return approximateLog10SumLog10(approximateLog10SumLog10(a, b), c);
    }

    final static double approximateLog10SumLog10(double small, double big) {
        // make sure small is really the smaller value
        if ( small > big ) {
            final double t = big;
            big = small;
            small = t;
        }

        if (small == Double.NEGATIVE_INFINITY || big == Double.NEGATIVE_INFINITY )
            return big;

        if (big >= small + MathUtils.MAX_JACOBIAN_TOLERANCE)
            return big;

        // OK, so |y-x| < tol: we use the following identity then:
        // we need to compute log10(10^x + 10^y)
        // By Jacobian logarithm identity, this is equal to
        // max(x,y) + log10(1+10^-abs(x-y))
        // we compute the second term as a table lookup
        // with integer quantization
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        //final int ind = (int)(((big-small)/JACOBIAN_LOG_TABLE_STEP)); // hard rounding
        int ind = (int)(Math.round((big-small)/MathUtils.JACOBIAN_LOG_TABLE_STEP)); // hard rounding

        //double z =Math.log10(1+Math.pow(10.0,-diff));
        //System.out.format("x: %f, y:%f, app: %f, true: %f ind:%d\n",x,y,t2,z,ind);
        return big + MathUtils.jacobianLogTable[ind];
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

        // todo = can't deal with optimal dynamic programming solution with multiallelic records
        if (SIMPLE_GREEDY_GENOTYPER || !vc.isBiallelic()) {
            sampleIndices.addAll(GLs.keySet());
            sampleIdx = GLs.size();
        }
        else {

            for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
                if ( !sample.getValue().hasLikelihoods() )
                    continue;

                double[] likelihoods = sample.getValue().getLikelihoods().getAsVector();

                if (MathUtils.sum(likelihoods) > SUM_GL_THRESH_NOCALL)     {
                    //System.out.print(sample.getKey()+":");
                    //for (int k=0; k < likelihoods.length; k++)
                    //   System.out.format("%4.2f ",likelihoods[k]);
                    //System.out.println();
                    // all likelihoods are essentially the same: skip this sample and will later on force no call.
                    //sampleIdx++;
                    continue;
                }

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
            // if all likelihoods are essentially the same: we want to force no-call. In this case, we skip this sample for now,
            // and will add no-call genotype to GL's in a second pass
            ArrayList<Allele> myAlleles = new ArrayList<Allele>();

            double qual = Double.NEGATIVE_INFINITY;
            double[] likelihoods = g.getLikelihoods().getAsVector();

            if (SIMPLE_GREEDY_GENOTYPER || !vc.isBiallelic()) {
                bestGTguess = Utils.findIndexOfMaxEntry(g.getLikelihoods().getAsVector());
            }
            else {
                int newIdx = tracebackArray[k][startIdx];;
                bestGTguess = startIdx - newIdx;
                startIdx = newIdx;
            }

            /*           System.out.format("Sample: %s GL:",sample);
                    for (int i=0; i < likelihoods.length; i++)
                        System.out.format("%1.4f, ",likelihoods[i]);
            */

            for (int i=0; i < likelihoods.length; i++) {
                if (i==bestGTguess)
                    continue;
                if (likelihoods[i] >= qual)
                    qual = likelihoods[i];
            }
            // qual contains now max(likelihoods[k]) for all k != bestGTguess
            qual = likelihoods[bestGTguess] - qual;

            // likelihoods are stored row-wise in lower triangular matrix. IE
            // for 2 alleles they have ordering AA,AB,BB
            // for 3 alleles they are ordered AA,AB,BB,AC,BC,CC
            // Get now alleles corresponding to best index
            int kk=0;
            boolean done = false;
            for (int j=0; j < vc.getNAlleles(); j++) {
                for (int i=0; i <= j; i++){
                    if (kk++ == bestGTguess) {
                        if (i==0)
                            myAlleles.add(vc.getReference());
                        else
                            myAlleles.add(vc.getAlternateAllele(i-1));

                        if (j==0)
                            myAlleles.add(vc.getReference());
                        else
                            myAlleles.add(vc.getAlternateAllele(j-1));
                        done = true;
                        break;
                    }

                }
                if (done)
                    break;
            }

            if (qual < 0) {
                // QUAL can be negative if the chosen genotype is not the most likely one individually.
                // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen on
                double[] normalized = MathUtils.normalizeFromLog10(likelihoods);
                double chosenGenotype = normalized[bestGTguess];
                qual = -1.0 * Math.log10(1.0 - chosenGenotype);
            }
            //System.out.println(myAlleles.toString());
            calls.put(sample, new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));

        }

        for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {

            if ( !sample.getValue().hasLikelihoods() )
                continue;
            Genotype g = GLs.get(sample.getKey());

            double[] likelihoods = sample.getValue().getLikelihoods().getAsVector();

            if (MathUtils.sum(likelihoods) <= SUM_GL_THRESH_NOCALL)
                continue; // regular likelihoods

            ArrayList<Allele> myAlleles = new ArrayList<Allele>();

            double qual = Genotype.NO_NEG_LOG_10PERROR;
            myAlleles.add(Allele.NO_CALL);
            myAlleles.add(Allele.NO_CALL);
            //System.out.println(myAlleles.toString());
            calls.put(sample.getKey(), new Genotype(sample.getKey(), myAlleles, qual, null, g.getAttributes(), false));
        }
        return calls;
    }

    // -------------------------------------------------------------------------------------
    //
    // Gold standard, but O(N^2), implementation.
    //
    // TODO -- remove me for clarity in this code
    //
    // -------------------------------------------------------------------------------------
    public int gdaN2GoldStandard(Map<String, Genotype> GLs,
                                 double[] log10AlleleFrequencyPriors,
                                 double[] log10AlleleFrequencyPosteriors, int idxAA, int idxAB, int idxBB) {
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
            double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];

            // special treatment for k=0: iteration reduces to:
            //YMatrix[j][0] = YMatrix[j-1][0]*genotypeLikelihoods[GenotypeType.AA.ordinal()];
            logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[idxAA];

            for (int k=1; k <= 2*j; k++ ) {

                //double num = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[GenotypeType.AA.ordinal()];
                double logNumerator[];
                logNumerator = new double[3];
                if (k < 2*j-1)
                    logNumerator[0] = MathUtils.log10Cache[2*j-k] + MathUtils.log10Cache[2*j-k-1] + logYMatrix[j-1][k] +
                            genotypeLikelihoods[idxAA];
                else
                    logNumerator[0] = Double.NEGATIVE_INFINITY;


                if (k < 2*j)
                    logNumerator[1] = MathUtils.log10Cache[2*k] + MathUtils.log10Cache[2*j-k]+ logYMatrix[j-1][k-1] +
                            genotypeLikelihoods[idxAB];
                else
                    logNumerator[1] = Double.NEGATIVE_INFINITY;

                if (k > 1)
                    logNumerator[2] = MathUtils.log10Cache[k] + MathUtils.log10Cache[k-1] + logYMatrix[j-1][k-2] +
                            genotypeLikelihoods[idxBB];
                else
                    logNumerator[2] = Double.NEGATIVE_INFINITY;

                double logNum = MathUtils.softMax(logNumerator);

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

}
