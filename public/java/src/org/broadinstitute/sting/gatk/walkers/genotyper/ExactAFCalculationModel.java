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
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {

    // private final static boolean DEBUG = false;

    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
    }

    public List<Allele> getLog10PNonRef(final VariantContext vc,
                                        final double[] log10AlleleFrequencyPriors,
                                        final AlleleFrequencyCalculationResult result) {

        GenotypesContext GLs = vc.getGenotypes();
        List<Allele> alleles = vc.getAlleles();

        final int myMaxAltAllelesToGenotype = CAP_MAX_ALTERNATE_ALLELES_FOR_INDELS && vc.getType().equals(VariantContext.Type.INDEL) ? 2 : MAX_ALTERNATE_ALLELES_TO_GENOTYPE;

        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > myMaxAltAllelesToGenotype ) {
            logger.warn("this tool is currently set to genotype at most " + myMaxAltAllelesToGenotype + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + (vc.getAlternateAlleles().size()) + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            alleles = new ArrayList<Allele>(myMaxAltAllelesToGenotype + 1);
            alleles.add(vc.getReference());
            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, myMaxAltAllelesToGenotype));
            GLs = VariantContextUtils.subsetDiploidAlleles(vc, alleles, false);
        }

        linearExactMultiAllelic(GLs, alleles.size() - 1, log10AlleleFrequencyPriors, result);

        return alleles;
    }


    private static final int PL_INDEX_OF_HOM_REF = 0;
    private static final List<Allele> chooseMostLikelyAlternateAlleles(VariantContext vc, int numAllelesToChoose) {
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ )
            likelihoodSums[i] = new LikelihoodSum(vc.getAlternateAllele(i));

        // based on the GLs, find the alternate alleles with the most probability; sum the GLs for the most likely genotype
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes());
        for ( final double[] likelihoods : GLs ) {
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PL_INDEX_OF_HOM_REF ) {
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                if ( alleles.alleleIndex1 != 0 )
                    likelihoodSums[alleles.alleleIndex1-1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                // don't double-count it
                if ( alleles.alleleIndex2 != 0 && alleles.alleleIndex2 != alleles.alleleIndex1 )
                    likelihoodSums[alleles.alleleIndex2-1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
            }
        }

        // sort them by probability mass and choose the best ones
        Collections.sort(Arrays.asList(likelihoodSums));
        final ArrayList<Allele> bestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( int i = 0; i < numAllelesToChoose; i++ )
            bestAlleles.add(likelihoodSums[i].allele);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( Allele allele : vc.getAlternateAlleles() ) {
            if ( bestAlleles.contains(allele) )
                orderedBestAlleles.add(allele);
        }
        
        return orderedBestAlleles;
    }


    // -------------------------------------------------------------------------------------
    //
    // Multi-allelic implementation.
    //
    // -------------------------------------------------------------------------------------

    public static void linearExactMultiAllelic(final GenotypesContext GLs,
                                               final int numAlternateAlleles,
                                               final double[] log10AlleleFrequencyPriors,
                                               final AlleleFrequencyCalculationResult result) {

        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final LinkedList<ExactACset> ACqueue = new LinkedList<ExactACset>();

        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<ExactACcounts, ExactACset>(numChr+1);

        // add AC=0 to the queue
        int[] zeroCounts = new int[numAlternateAlleles];
        ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.ACcounts, zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        double maxLog10L = Double.NEGATIVE_INFINITY;
        while ( !ACqueue.isEmpty() ) {
            // compute log10Likelihoods
            final ExactACset set = ACqueue.remove();
            final double log10LofKs = calculateAlleleCountConformation(set, genotypeLikelihoods, maxLog10L, numChr, ACqueue, indexesToACset, log10AlleleFrequencyPriors, result);

            // adjust max likelihood seen if needed
            maxLog10L = Math.max(maxLog10L, log10LofKs);

            // clean up memory
            indexesToACset.remove(set.ACcounts);
            //if ( DEBUG )
            //    System.out.printf(" *** removing used set=%s%n", set.ACcounts);
        }
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;
        
        public DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }

    private static double calculateAlleleCountConformation(final ExactACset set,
                                                           final ArrayList<double[]> genotypeLikelihoods,
                                                           final double maxLog10L,
                                                           final int numChr,
                                                           final LinkedList<ExactACset> ACqueue,
                                                           final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                           final double[] log10AlleleFrequencyPriors,
                                                           final AlleleFrequencyCalculationResult result) {

        //if ( DEBUG )
        //    System.out.printf(" *** computing LofK for set=%s%n", set.ACcounts);

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, log10AlleleFrequencyPriors, result);

        final double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
            //if ( DEBUG )
            //    System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.ACcounts, log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;

        final int numAltAlleles = set.ACcounts.getCounts().length;

        // add conformations for the k+1 case
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.ACcounts.getCounts().clone();
            ACcountsClone[allele]++;
            // to get to this conformation, a sample would need to be AB (remember that ref=0)
            final int PLindex = GenotypeLikelihoods.calculatePLindex(0, allele+1);
            updateACset(ACcountsClone, numChr, set, PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final ArrayList<DependentSet> differentAlleles = new ArrayList<DependentSet>(numAltAlleles * numAltAlleles);
            final ArrayList<DependentSet> sameAlleles = new ArrayList<DependentSet>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.ACcounts.getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    // to get to this conformation, a sample would need to be BB or BC (remember that ref=0, so add one to the index)
                    final int PLindex = GenotypeLikelihoods.calculatePLindex(allele_i+1, allele_j+1);
                    if ( allele_i == allele_j )
                        sameAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    else
                        differentAlleles.add(new DependentSet(ACcountsClone, PLindex));
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( DependentSet dependent : differentAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            for ( DependentSet dependent : sameAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        return log10LofK;
    }

    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also pushes its value to the given callingSetIndex.
    private static void updateACset(final int[] newSetCounts,
                                    final int numChr,
                                    final ExactACset dependentSet,
                                    final int PLsetIndex,
                                    final Queue<ExactACset> ACqueue,
                                    final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                    final ArrayList<double[]> genotypeLikelihoods) {
        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // push data from the dependency to the new set
        //if ( DEBUG )
        //    System.out.println(" *** pushing data from " + index + " to " + dependencySet.ACcounts);
        pushData(indexesToACset.get(index), dependentSet, PLsetIndex, genotypeLikelihoods);
    }

    private static void computeLofK(final ExactACset set,
                                    final ArrayList<double[]> genotypeLikelihoods,
                                    final double[] log10AlleleFrequencyPriors,
                                    final AlleleFrequencyCalculationResult result) {

        set.log10Likelihoods[0] = 0.0; // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for ( int j = 1; j < set.log10Likelihoods.length; j++ )
                set.log10Likelihoods[j] = set.log10Likelihoods[j-1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];

            final double log10Lof0 = set.log10Likelihoods[set.log10Likelihoods.length-1];
            result.setLog10LikelihoodOfAFzero(log10Lof0);
            result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return;
        }

        // if we got here, then k > 0 for at least one k.
        // the non-AA possible conformations were already dealt with by pushes from dependent sets;
        // now deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
        for ( int j = 1; j < set.log10Likelihoods.length; j++ ) {

            if ( totalK < 2*j-1 ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue = MathUtils.log10Cache[2*j-totalK] + MathUtils.log10Cache[2*j-totalK-1] + set.log10Likelihoods[j-1] + gl[HOM_REF_INDEX];
                set.log10Likelihoods[j] = MathUtils.approximateLog10SumLog10(set.log10Likelihoods[j], conformationValue);
            }

            final double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];
            set.log10Likelihoods[j] = set.log10Likelihoods[j] - logDenominator;
        }

        double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // update the MLE if necessary
        result.updateMLEifNeeded(log10LofK, set.ACcounts.counts);

        // apply the priors over each alternate allele
        for ( final int ACcount : set.ACcounts.getCounts() ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        result.updateMAPifNeeded(log10LofK, set.ACcounts.counts);
    }

    private static void pushData(final ExactACset targetSet,
                                 final ExactACset dependentSet,
                                 final int PLsetIndex,
                                 final ArrayList<double[]> genotypeLikelihoods) {
        final int totalK = targetSet.getACsum();

        for ( int j = 1; j < targetSet.log10Likelihoods.length; j++ ) {

            if ( totalK <= 2*j ) { // skip impossible conformations
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue =
                        determineCoefficient(PLsetIndex, j, targetSet.ACcounts.getCounts(), totalK) + dependentSet.log10Likelihoods[j-1] + gl[PLsetIndex];
                targetSet.log10Likelihoods[j] = MathUtils.approximateLog10SumLog10(targetSet.log10Likelihoods[j], conformationValue);
            }
        }        
    }

    private static double determineCoefficient(int PLindex, final int j, final int[] ACcounts, final int totalK) {

        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        // find the 2 alleles that are represented by this PL index
        GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);

        // *** note that throughout this method we subtract one from the alleleIndex because ACcounts ***
        // *** doesn't consider the reference allele whereas the GenotypeLikelihoods PL cache does.   ***

        // the AX het case
        if ( alleles.alleleIndex1 == 0 )
            return MathUtils.log10Cache[2*ACcounts[alleles.alleleIndex2-1]] + MathUtils.log10Cache[2*j-totalK];

        final int k_i = ACcounts[alleles.alleleIndex1-1];

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles.alleleIndex1 == alleles.alleleIndex2 ) {
            coeff = MathUtils.log10Cache[k_i] + MathUtils.log10Cache[k_i - 1];
        }
        // the het non-ref case (e.g. BC, BD, CD)
        else {
            final int k_j = ACcounts[alleles.alleleIndex2-1];
            coeff = MathUtils.log10Cache[2] + MathUtils.log10Cache[k_i] + MathUtils.log10Cache[k_j];
        }

        return coeff;
    }

    public GenotypesContext subsetAlleles(final VariantContext vc,
                                                      final List<Allele> allelesToUse,
                                                      final boolean assignGenotypes,
                                                      final int ploidy) {
        return VariantContextUtils.subsetDiploidAlleles(vc, allelesToUse, assignGenotypes);
    }

    // -------------------------------------------------------------------------------------
    //
    // Deprecated bi-allelic ~O(N) implementation.  Kept here for posterity.
    //
    // -------------------------------------------------------------------------------------

    /**
     * A simple data structure that holds the current, prev, and prev->prev likelihoods vectors
     * for the exact model calculation
     */
/*
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
                           double[][] log10AlleleFrequencyLikelihoods,
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
            log10AlleleFrequencyLikelihoods[0][k] = log10LofK;
            log10AlleleFrequencyPosteriors[0][k] = log10LofK + log10AlleleFrequencyPriors[k];

            // can we abort early?
            lastK = k;
            maxLog10L = Math.max(maxLog10L, log10LofK);
            if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
                //if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
                done = true;
            }

            logY.rotate();
        }

        return lastK;
    }

    final static double approximateLog10SumLog10(double a, double b, double c) {
        return approximateLog10SumLog10(approximateLog10SumLog10(a, b), c);
    }
*/

}
