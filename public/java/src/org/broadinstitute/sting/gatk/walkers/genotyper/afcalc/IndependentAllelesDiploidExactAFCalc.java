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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Computes the conditional bi-allelic exact results
 *
 * Suppose vc contains 2 alt allele: A* with C and T.  This function first computes:
 *
 * (1) P(D | AF_c > 0 && AF_t == *) [i.e., T can be anything]
 *
 * it then computes the conditional probability on AF_c == 0:
 *
 * (2) P(D | AF_t > 0 && AF_c == 0)
 *
 * Thinking about this visually, we have the following likelihood matrix where each cell is
 * the P(D | AF_c == i && AF_t == j):
 *
 *     0 AF_c > 0
 *    -----------------
 * 0  |  |
 *    |--|-------------
 * a  |  |
 * f  |  |
 * _  |  |
 * t  |  |
 * >  |  |
 * 0  |  |
 *
 * What we really want to know how
 *
 * (3) P(D | AF_c == 0 & AF_t == 0)
 *
 * compares with
 *
 * (4) P(D | AF_c > 0 || AF_t > 0)
 *
 * This is effectively asking for the value in the upper left vs. the sum of all cells.
 *
 * This class implements the conditional likelihoods summation for any number of alt
 * alleles, where each alt allele has its EXACT probability of segregating calculated by
 * reducing each alt B into the case XB and computing P(D | AF_b > 0 ) as follows:
 *
 * Suppose we have for a A/B/C site the following GLs:
 *
 * AA AB BB AC BC CC
 *
 * and we want to get the bi-allelic GLs for X/B, where X is everything not B
 *
 * XX = AA + AC + CC (since X = A or C)
 * XB = AB + BC
 * BB = BB
 *
 * After each allele has its probability calculated we compute the joint posterior
 * as P(D | AF_* == 0) = prod_i P (D | AF_i == 0), after applying the theta^i
 * prior for the ith least likely allele.
 */
 public class IndependentAllelesDiploidExactAFCalc extends DiploidExactAFCalc {
    /**
     * The min. confidence of an allele to be included in the joint posterior.
     */
    private final static double MIN_LOG10_CONFIDENCE_TO_INCLUDE_ALLELE_IN_POSTERIOR = Math.log10(1e-10);

    private final static int[] BIALLELIC_NON_INFORMATIVE_PLS = new int[]{0,0,0};
    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Sorts AFCalcResults by their posteriors of AF > 0, so the
     */
    private final static class CompareAFCalcResultsByPNonRef implements Comparator<AFCalcResult> {
        @Override
        public int compare(AFCalcResult o1, AFCalcResult o2) {
            return -1 * Double.compare(o1.getLog10PosteriorOfAFGT0(), o2.getLog10PosteriorOfAFGT0());
        }
    }

    private final static CompareAFCalcResultsByPNonRef compareAFCalcResultsByPNonRef = new CompareAFCalcResultsByPNonRef();

    /**
     * The AFCalc model we are using to do the bi-allelic computation
     */
    final AFCalc biAlleleExactModel;

    protected IndependentAllelesDiploidExactAFCalc(int nSamples, int maxAltAlleles, final int ploidy) {
        super(nSamples, maxAltAlleles, ploidy);
        biAlleleExactModel = new ReferenceDiploidExactAFCalc(nSamples, 1, ploidy);
    }

    /**
     * Trivial subclass that helps with debugging by keeping track of the supporting information for this joint call
     */
    private static class MyAFCalcResult extends AFCalcResult {
        /**
         * List of the supporting bi-allelic AFCalcResults that went into making this multi-allelic joint call
         */
        final List<AFCalcResult> supporting;

        private MyAFCalcResult(int[] alleleCountsOfMLE, int nEvaluations, List<Allele> allelesUsedInGenotyping, double[] log10LikelihoodsOfAC, double[] log10PriorsOfAC, Map<Allele, Double> log10pRefByAllele, List<AFCalcResult> supporting) {
            super(alleleCountsOfMLE, nEvaluations, allelesUsedInGenotyping, log10LikelihoodsOfAC, log10PriorsOfAC, log10pRefByAllele);
            this.supporting = supporting;
        }
    }

    @Override
    public AFCalcResult computeLog10PNonRef(final VariantContext vc,
                                            final double[] log10AlleleFrequencyPriors) {
        final List<AFCalcResult> independentResultTrackers = computeAlleleIndependentExact(vc, log10AlleleFrequencyPriors);
        final List<AFCalcResult> withMultiAllelicPriors = applyMultiAllelicPriors(independentResultTrackers);
        return combineIndependentPNonRefs(vc, withMultiAllelicPriors);
    }


    /**
     * Compute the conditional exact AFCalcResult for each allele in vc independently, returning
     * the result of each, in order of the alt alleles in VC
     *
     * @param vc the VariantContext we want to analyze
     * @param log10AlleleFrequencyPriors the priors
     * @return a list of the AFCalcResults for each bi-allelic sub context of vc
     */
    @Requires({"vc != null", "log10AlleleFrequencyPriors != null"})
    @Ensures("goodIndependentResult(vc, result)")
    protected final List<AFCalcResult> computeAlleleIndependentExact(final VariantContext vc,
                                                                     final double[] log10AlleleFrequencyPriors) {
        final List<AFCalcResult> results = new LinkedList<AFCalcResult>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc) ) {
            final AFCalcResult resultTracker = biAlleleExactModel.getLog10PNonRef(subvc, log10AlleleFrequencyPriors);
            results.add(resultTracker);
        }

        return results;
    }

    /**
     * Helper function to ensure that the computeAlleleIndependentExact is returning reasonable results
     */
    private static boolean goodIndependentResult(final VariantContext vc, final List<AFCalcResult> results) {
        if ( results.size() != vc.getNAlleles() - 1) return false;
        for ( int i = 0; i < results.size(); i++ ) {
            if ( results.get(i).getAllelesUsedInGenotyping().size() != 2 )
                return false;
            if ( ! results.get(i).getAllelesUsedInGenotyping().contains(vc.getAlternateAllele(i)) )
                return false;
        }

        return true;
    }

    /**
     * Returns the bi-allelic variant context for each alt allele in vc with bi-allelic likelihoods, in order
     *
     * @param vc the variant context to split.  Must have n.alt.alleles > 1
     * @return a bi-allelic variant context for each alt allele in vc
     */
    @Requires({"vc != null", "vc.getNAlleles() > 1"})
    @Ensures("result.size() == vc.getNAlleles() - 1")
    protected final List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;
        final List<VariantContext> vcs = new LinkedList<VariantContext>();

        for ( int altI = 0; altI < nAltAlleles; altI++ ) {
            vcs.add(biallelicCombinedGLs(vc, altI + 1));
        }

        return vcs;
    }

    /**
     * Create a single bi-allelic variant context from rootVC with alt allele with index altAlleleIndex
     *
     * @param rootVC the root (potentially multi-allelic) variant context
     * @param altAlleleIndex index of the alt allele, from 0 == first alt allele
     * @return a bi-allelic variant context based on rootVC
     */
    @Requires({"rootVC.getNAlleles() > 1", "altAlleleIndex < rootVC.getNAlleles()"})
    @Ensures({"result.isBiallelic()"})
    protected final VariantContext biallelicCombinedGLs(final VariantContext rootVC, final int altAlleleIndex) {
        if ( rootVC.isBiallelic() ) {
            return rootVC;
        } else {
            final int nAlts = rootVC.getNAlleles() - 1;
            final List<Genotype> biallelicGenotypes = new ArrayList<Genotype>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLs(g, altAlleleIndex, nAlts));

            final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
            final Allele altAllele = rootVC.getAlternateAllele(altAlleleIndex - 1);
            vcb.alleles(Arrays.asList(rootVC.getReference(), altAllele));
            vcb.genotypes(biallelicGenotypes);
            return vcb.make();
        }
    }

    /**
     * Returns a new Genotype with the PLs of the multi-allelic original reduced to a bi-allelic case
     *
     * This is handled in the following way:
     *
     * Suppose we have for a A/B/C site the following GLs:
     *
     * AA AB BB AC BC CC
     *
     * and we want to get the bi-allelic GLs for X/B, where X is everything not B
     *
     * XX = AA + AC + CC (since X = A or C)
     * XB = AB + BC
     * BB = BB
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @Requires({"original.hasLikelihoods()"}) // TODO -- add ploidy == 2 test "original.getPLs() == null || original.getPLs().length == 3"})
    @Ensures({"result.hasLikelihoods()", "result.getPL().length == 3"})
    protected Genotype combineGLs(final Genotype original, final int altIndex, final int nAlts ) {
        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(BIALLELIC_NON_INFORMATIVE_PLS).alleles(BIALLELIC_NOCALL).make();

        if ( altIndex < 1 || altIndex > nAlts ) throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);

        final double[] normalizedPr = MathUtils.normalizeFromLog10(GenotypeLikelihoods.fromPLs(original.getPL()).getAsVector());
        final double[] biAllelicPr = new double[3];

        for ( int index = 0; index < normalizedPr.length; index++ ) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(index);

            if ( pair.alleleIndex1 == altIndex ) {
                if ( pair.alleleIndex2 == altIndex )
                    // hom-alt case
                    biAllelicPr[2] = normalizedPr[index];
                else
                    // het-alt case
                    biAllelicPr[1] += normalizedPr[index];
            } else {
                if ( pair.alleleIndex2 == altIndex )
                    // het-alt case
                    biAllelicPr[1] += normalizedPr[index];
                else
                    // hom-non-alt
                    biAllelicPr[0] += normalizedPr[index];
            }
        }

        final double[] GLs = new double[3];
        for ( int i = 0; i < GLs.length; i++ ) GLs[i] = Math.log10(biAllelicPr[i]);

        return new GenotypeBuilder(original).PL(GLs).alleles(BIALLELIC_NOCALL).make();
    }

    protected final List<AFCalcResult> applyMultiAllelicPriors(final List<AFCalcResult> conditionalPNonRefResults) {
        final ArrayList<AFCalcResult> sorted = new ArrayList<AFCalcResult>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, compareAFCalcResultsByPNonRef);

        double lastPosteriorGt0 = sorted.get(0).getLog10PosteriorOfAFGT0();
        final double log10SingleAllelePriorOfAFGt0 = conditionalPNonRefResults.get(0).getLog10PriorOfAFGT0();

        for ( int i = 0; i < sorted.size(); i++ ) {
            if ( sorted.get(i).getLog10PosteriorOfAFGT0() > lastPosteriorGt0 )
                throw new IllegalStateException("pNonRefResults not sorted: lastPosteriorGt0 " + lastPosteriorGt0 + " but current is " + sorted.get(i).getLog10PosteriorOfAFGT0());

                final double log10PriorAFGt0 = (i + 1) * log10SingleAllelePriorOfAFGt0;
            final double log10PriorAFEq0 = Math.log10(1 - Math.pow(10, log10PriorAFGt0));
            final double[] thetaTONPriors = new double[] { log10PriorAFEq0, log10PriorAFGt0 };

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            sorted.set(i, sorted.get(i).withNewPriors(MathUtils.normalizeFromLog10(thetaTONPriors, true)));
        }

        return sorted;
    }


    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Given n independent calculations for each of n alternate alleles create a single
     * combined AFCalcResult with:
     *
     * priors for AF == 0 equal to theta^N for the nth least likely allele
     * posteriors that reflect the combined chance that any alleles are segregating and corresponding
     * likelihoods
     * combined MLEs in the order of the alt alleles in vc
     *
     * @param sortedResultsWithThetaNPriors the pNonRef result for each allele independently
     */
    protected AFCalcResult combineIndependentPNonRefs(final VariantContext vc,
                                                      final List<AFCalcResult> sortedResultsWithThetaNPriors) {
        int nEvaluations = 0;
        final int nAltAlleles = sortedResultsWithThetaNPriors.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final double[] log10PriorsOfAC = new double[2];
        final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

        // the sum of the log10 posteriors for AF == 0 and AF > 0 to determine joint probs
        double log10PosteriorOfACEq0Sum = 0.0;
        double log10PosteriorOfACGt0Sum = 0.0;

        boolean anyPoly = false;
        for ( final AFCalcResult sortedResultWithThetaNPriors : sortedResultsWithThetaNPriors ) {
            final Allele altAllele = sortedResultWithThetaNPriors.getAllelesUsedInGenotyping().get(1);
            final int altI = vc.getAlleles().indexOf(altAllele) - 1;

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = sortedResultWithThetaNPriors.getAlleleCountAtMLE(altAllele);

            // the AF > 0 case requires us to store the normalized likelihood for later summation
            if ( sortedResultWithThetaNPriors.getLog10PosteriorOfAFGT0() > MIN_LOG10_CONFIDENCE_TO_INCLUDE_ALLELE_IN_POSTERIOR ) {
                anyPoly = true;
                log10PosteriorOfACEq0Sum += sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0();
                log10PriorsOfAC[0] += sortedResultWithThetaNPriors.getLog10PriorOfAFEq0();
                log10PriorsOfAC[1] += sortedResultWithThetaNPriors.getLog10PriorOfAFGT0();
            }

            log10PosteriorOfACGt0Sum += sortedResultWithThetaNPriors.getLog10PosteriorOfAFGT0();

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            log10pRefByAllele.put(altAllele, sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0());

            // trivial -- update the number of evaluations
            nEvaluations += sortedResultWithThetaNPriors.nEvaluations;
        }

        // If no alleles were polymorphic, make sure we have the proper priors (the defaults) for likelihood calculation
        if ( ! anyPoly ) {
            log10PriorsOfAC[0] = sortedResultsWithThetaNPriors.get(0).getLog10PriorOfAFEq0();
            log10PriorsOfAC[1] = sortedResultsWithThetaNPriors.get(0).getLog10PriorOfAFGT0();
        }

        // In principle, if B_p = x and C_p = y are the probabilities of being poly for alleles B and C,
        // the probability of being poly is (1 - B_p) * (1 - C_p) = (1 - x) * (1 - y).  We want to estimate confidently
        // log10((1 - x) * (1 - y)) which is log10(1 - x) + log10(1 - y).  This sum is log10PosteriorOfACEq0
        //
        // note we need to handle the case where the posterior of AF == 0 is 0.0, in which case we
        // use the summed log10PosteriorOfACGt0Sum directly.  This happens in cases where
        //   AF > 0 : 0.0 and AF == 0 : -16, and if you use the inverse calculation you get 0.0 and MathUtils.LOG10_P_OF_ZERO
        final double log10PosteriorOfACGt0;
        if ( log10PosteriorOfACEq0Sum == 0.0 )
            log10PosteriorOfACGt0 = log10PosteriorOfACGt0Sum;
        else
            log10PosteriorOfACGt0 = Math.max(Math.log10(1 - Math.pow(10, log10PosteriorOfACEq0Sum)), MathUtils.LOG10_P_OF_ZERO);

        final double[] log10LikelihoodsOfAC = new double[] {
                // L + prior = posterior => L = poster - prior
                log10PosteriorOfACEq0Sum - log10PriorsOfAC[0],
                log10PosteriorOfACGt0 - log10PriorsOfAC[1]
        };

        return new MyAFCalcResult(alleleCountsOfMLE, nEvaluations, vc.getAlleles(),
                // necessary to ensure all values < 0
                MathUtils.normalizeFromLog10(log10LikelihoodsOfAC, true),
                // priors incorporate multiple alt alleles, must be normalized
                MathUtils.normalizeFromLog10(log10PriorsOfAC, true),
                log10pRefByAllele, sortedResultsWithThetaNPriors);
    }
}
