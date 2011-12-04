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
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {
    //
    // code for testing purposes
    //
    private final static boolean DEBUG = false;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
    private final static double SUM_GL_THRESH_NOCALL = -0.001; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    private static final boolean SIMPLE_GREEDY_GENOTYPER = false;
    private static final List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    private final boolean USE_MULTI_ALLELIC_CALCULATION;


    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        USE_MULTI_ALLELIC_CALCULATION = UAC.MULTI_ALLELIC;
    }

    public void getLog10PNonRef(GenotypesContext GLs, List<Allele> alleles,
                                double[] log10AlleleFrequencyPriors,
                                double[][] log10AlleleFrequencyPosteriors) {
        final int numAlleles = alleles.size();

        if ( USE_MULTI_ALLELIC_CALCULATION )
            linearExactMultiAllelic(GLs, numAlleles - 1, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors, false);
        else
            linearExact(GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);
    }

    private static final ArrayList<double[]> getGLs(GenotypesContext GLs) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>();

        genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
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

    public int linearExact(GenotypesContext GLs,
                           double[] log10AlleleFrequencyPriors,
                           double[][] log10AlleleFrequencyPosteriors) {
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
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods.get(j)[0];
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
                        aa = MathUtils.log10Cache[2*j-k] + MathUtils.log10Cache[2*j-k-1] + kMinus0[j-1] + gl[0];

                    if (k < 2*j)
                        ab = MathUtils.log10Cache[2*k] + MathUtils.log10Cache[2*j-k]+ kMinus1[j-1] + gl[1];

                    double log10Max;
                    if (k > 1) {
                        final double bb = MathUtils.log10Cache[k] + MathUtils.log10Cache[k-1] + kMinus2[j-1] + gl[2];
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
            log10AlleleFrequencyPosteriors[0][k] = log10LofK + log10AlleleFrequencyPriors[k];

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

    final static double approximateLog10SumLog10(double[] vals) {
        if ( vals.length < 2 )
            throw new ReviewedStingException("Passing array with fewer than 2 values when computing approximateLog10SumLog10");

        double approx = approximateLog10SumLog10(vals[0], vals[1]);
        for ( int i = 2; i < vals.length; i++ )
            approx = approximateLog10SumLog10(approx, vals[i]);
        return approx;
    }

    final static double approximateLog10SumLog10(double a, double b, double c) {
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


    // -------------------------------------------------------------------------------------
    //
    // Multi-allelic implementation.
    //
    // -------------------------------------------------------------------------------------

    private static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first
    private static final int AC_ZERO_INDEX = 0;  // ExactACset index for k=0 over all k

    // This class represents a column in the Exact AC calculation matrix
    private static final class ExactACset {

        // the counts of the various alternate alleles which this column represents
        final int[] ACcounts;

        // the column of the matrix
        final double[] log10Likelihoods;

        // mapping of column index for those columns upon which this one depends to the index into the PLs which is used as the transition to this column;
        // for example, in the biallelic case, the transition from k=0 to k=1 would be AB while the transition to k=2 would be BB.
        final HashMap<Integer, Integer> ACsetIndexToPLIndex = new HashMap<Integer, Integer>();

        // to minimize memory consumption, we know we can delete any sets in this list because no further sets will depend on them
        final ArrayList<Integer> dependentACsetsToDelete = new ArrayList<Integer>();

        // index used to represent this set in the global hashmap: (numSamples^0 * allele_1) + (numSamples^1 * allele_2) + (numSamples^2 * allele_3) + ...
        private int index = -1;

        public ExactACset(int size, int[] ACcounts) {
            this.ACcounts = ACcounts;
            log10Likelihoods = new double[size];
        }

        public int getIndex() {
            if ( index == -1 )
                index = generateIndex(ACcounts, log10Likelihoods.length);
            return index;
        }

        public static int generateIndex(int[] ACcounts, int multiplier) {
            int index = 0;
            for ( int i = 0; i < ACcounts.length; i++ )
                index += Math.pow(multiplier, i) * ACcounts[i];
            return index;
        }

        // sum of all the non-reference alleles
        public int getACsum() {
            int sum = 0;
            for ( int count : ACcounts )
                sum += count;
            return sum;
        }
    }

    public static void linearExactMultiAllelic(GenotypesContext GLs,
                                               int numAlternateAlleles,
                                               double[] log10AlleleFrequencyPriors,
                                               double[][] log10AlleleFrequencyPosteriors,
                                               boolean preserveData) {

        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final Queue<ExactACset> ACqueue = new LinkedList<ExactACset>();

        // mapping of ExactACset indexes to the objects
        final HashMap<Integer, ExactACset> indexesToACset = new HashMap<Integer, ExactACset>(numChr+1);

        // add AC=0 to the queue
        int[] zeroCounts = new int[numAlternateAlleles];
        ExactACset zeroSet = new ExactACset(numSamples+1, zeroCounts);
        ACqueue.add(zeroSet);
        indexesToACset.put(0, zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        double maxLog10L = Double.NEGATIVE_INFINITY;
        while ( !ACqueue.isEmpty() ) {
            // compute log10Likelihoods
            final ExactACset set = ACqueue.remove();
            final double log10LofKs = calculateAlleleCountConformation(set, genotypeLikelihoods, maxLog10L, numChr, preserveData, ACqueue, indexesToACset, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);

            // adjust max likelihood seen if needed
            maxLog10L = Math.max(maxLog10L, log10LofKs);
        }
    }

    private static double calculateAlleleCountConformation(final ExactACset set,
                                                           final ArrayList<double[]> genotypeLikelihoods,
                                                           final double maxLog10L,
                                                           final int numChr,
                                                           final boolean preserveData,
                                                           final Queue<ExactACset> ACqueue,
                                                           final HashMap<Integer, ExactACset> indexesToACset,
                                                           final double[] log10AlleleFrequencyPriors,
                                                           final double[][] log10AlleleFrequencyPosteriors) {

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, indexesToACset, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors);

        // clean up memory
        if ( !preserveData ) {
            for ( int index : set.dependentACsetsToDelete )
                indexesToACset.put(index, null);
        }

        final double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
            if ( DEBUG ) System.out.printf(" *** breaking early ks=%d log10L=%.2f maxLog10L=%.2f%n", set.index, log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;

        ExactACset lastSet = null; // keep track of the last set placed in the queue so that we can tell it to clean us up when done processing
        int numAltAlleles = set.ACcounts.length;

        // genotype likelihoods are a linear vector that can be thought of as a row-wise upper triangular matrix of log10Likelihoods.
        // so e.g. with 2 alt alleles the likelihoods are AA,AB,AC,BB,BC,CC and with 3 alt alleles they are AA,AB,AC,AD,BB,BC,BD,CC,CD,DD.

        // add conformations for the k+1 case
        int PLindex = 0;
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            int[] ACcountsClone = set.ACcounts.clone();
            ACcountsClone[allele]++;
            lastSet = updateACset(ACcountsClone, numChr, set.getIndex(), ++PLindex, ACqueue, indexesToACset);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    int[] ACcountsClone = set.ACcounts.clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;
                    lastSet = updateACset(ACcountsClone, numChr,set.getIndex(), ++PLindex , ACqueue, indexesToACset);
                }
            }
        }

        if ( lastSet == null )
            throw new ReviewedStingException("No new AC sets were added or updated but the AC still hasn't reached 2N");
        lastSet.dependentACsetsToDelete.add(set.index);

        return log10LofK;
    }

    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also adds it as a dependency to the given callingSetIndex.
    private static ExactACset updateACset(int[] ACcounts,
                                          int numChr,
                                          final int callingSetIndex,
                                          final int PLsetIndex,
                                          final Queue<ExactACset> ACqueue,
                                          final HashMap<Integer, ExactACset> indexesToACset) {
        final int index = ExactACset.generateIndex(ACcounts, numChr+1);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset set = new ExactACset(numChr/2 +1, ACcounts);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // add the given dependency to the set
        ExactACset set = indexesToACset.get(index);
        set.ACsetIndexToPLIndex.put(callingSetIndex, PLsetIndex);
        return set;
    }

    private static void computeLofK(ExactACset set,
                                    ArrayList<double[]> genotypeLikelihoods,
                                    final HashMap<Integer, ExactACset> indexesToACset,
                                    double[] log10AlleleFrequencyPriors,
                                    double[][] log10AlleleFrequencyPosteriors) {

        set.log10Likelihoods[0] = 0.0; // the zero case
        int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( set.getIndex() == AC_ZERO_INDEX ) {
            for ( int j = 1; j < set.log10Likelihoods.length; j++ )
                set.log10Likelihoods[j] = set.log10Likelihoods[j-1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];
        }
        // k > 0 for at least one k
        else {
            // all possible likelihoods for a given cell from which to choose the max
            final int numPaths = set.ACsetIndexToPLIndex.size() + 1;
            final double[] log10ConformationLikelihoods = new double[numPaths];

            for ( int j = 1; j < set.log10Likelihoods.length; j++ ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];

                // initialize
                for ( int i = 0; i < numPaths; i++ )
                    log10ConformationLikelihoods[i] = Double.NEGATIVE_INFINITY;

                // deal with the AA case first
                if ( totalK < 2*j-1 )
                    log10ConformationLikelihoods[0] = MathUtils.log10Cache[2*j-totalK] + MathUtils.log10Cache[2*j-totalK-1] + set.log10Likelihoods[j-1] + gl[HOM_REF_INDEX];

                // deal with the other possible conformations now
                if ( totalK <= 2*j ) { // skip impossible conformations
                    int conformationIndex = 1;
                    for ( Map.Entry<Integer, Integer> mapping : set.ACsetIndexToPLIndex.entrySet() )
                        log10ConformationLikelihoods[conformationIndex++] =
                                determineCoefficient(mapping.getValue(), j, totalK) + indexesToACset.get(mapping.getKey()).log10Likelihoods[j-1] + gl[mapping.getValue()];
                }

                double log10Max = approximateLog10SumLog10(log10ConformationLikelihoods);

                // finally, update the L(j,k) value
                set.log10Likelihoods[j] = log10Max - logDenominator;
            }
        }

        // update the posteriors vector
        final double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // determine the power of theta to use
        int nonRefAlleles = 0;
        for ( int i = 0; i < set.ACcounts.length; i++ ) {
            if ( set.ACcounts[i] > 0 )
                nonRefAlleles++;
        }
        if ( nonRefAlleles == 0 ) // for k=0 we still want to use a power of 1
            nonRefAlleles++;

        // update the posteriors vector which is a collapsed view of each of the various ACs
        for ( int i = 0; i < set.ACcounts.length; i++ ) {
            // TODO -- double check the math and then cache these values for efficiency
            double prior = Math.pow(log10AlleleFrequencyPriors[totalK], nonRefAlleles);
            log10AlleleFrequencyPosteriors[i][set.ACcounts[i]] = approximateLog10SumLog10(log10AlleleFrequencyPosteriors[i][set.ACcounts[i]], log10LofK + prior);
        }
    }

    private static double determineCoefficient(int PLindex, int j, int totalK) {

        // TODO -- the math here needs to be fixed and checked; hard-coding in the biallelic case
        //AA,AB,AC,AD,BB,BC,BD,CC,CD,DD.

        double coeff;
        if ( PLindex == 1 )
            coeff = MathUtils.log10Cache[2*totalK] + MathUtils.log10Cache[2*j-totalK];
        else
            coeff = MathUtils.log10Cache[totalK] + MathUtils.log10Cache[totalK-1];
        return coeff;
    }

    /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    public GenotypesContext assignGenotypes(VariantContext vc,
                                            double[][] log10AlleleFrequencyPosteriors,
                                            int AFofMaxLikelihood) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        GenotypesContext GLs = vc.getGenotypes();
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
            sampleIndices.addAll(GLs.getSampleNamesOrderedByName());
            sampleIdx = GLs.size();
        }
        else {

            for ( final Genotype genotype : GLs.iterateInSampleNameOrder() ) {
                if ( !genotype.hasLikelihoods() )
                    continue;

                double[] likelihoods = genotype.getLikelihoods().getAsVector();

                if (MathUtils.sum(likelihoods) > SUM_GL_THRESH_NOCALL)     {
                    //System.out.print(sample.getKey()+":");
                    //for (int k=0; k < likelihoods.length; k++)
                    //   System.out.format("%4.2f ",likelihoods[k]);
                    //System.out.println();
                    // all likelihoods are essentially the same: skip this sample and will later on force no call.
                    //sampleIdx++;
                    continue;
                }

                sampleIndices.add(genotype.getSampleName());

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

        GenotypesContext calls = GenotypesContext.create();

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

            double[] likelihoods = g.getLikelihoods().getAsVector();

            if (SIMPLE_GREEDY_GENOTYPER || !vc.isBiallelic()) {
                bestGTguess = Utils.findIndexOfMaxEntry(likelihoods);
            }
            else {
                int newIdx = tracebackArray[k][startIdx];;
                bestGTguess = startIdx - newIdx;
                startIdx = newIdx;
            }

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

            final double qual = GenotypeLikelihoods.getQualFromLikelihoods(bestGTguess, likelihoods);
            //System.out.println(myAlleles.toString());
            calls.add(new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));
        }

        for ( final Genotype genotype : GLs.iterateInSampleNameOrder() ) {
            if ( !genotype.hasLikelihoods() )
                continue;

            final Genotype g = GLs.get(genotype.getSampleName());
            final double[] likelihoods = genotype.getLikelihoods().getAsVector();

            if (MathUtils.sum(likelihoods) <= SUM_GL_THRESH_NOCALL)
                continue; // regular likelihoods

            final double qual = Genotype.NO_LOG10_PERROR;
            calls.replace(new Genotype(g.getSampleName(), NO_CALL_ALLELES, qual, null, g.getAttributes(), false));
        }

        return calls;
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
