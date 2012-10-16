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
    private final static double MIN_LOG10_CONFIDENCE_TO_INCLUDE_ALLELE_IN_POSTERIOR = Math.log10(1e-20);

    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Sorts AFCalcResults by their posteriors of AF > 0, so the
     */
    private final static class CompareAFCalcResultsByPNonRef implements Comparator<AFCalcResult> {
        @Override
        public int compare(AFCalcResult o1, AFCalcResult o2) {
            return Double.compare(o1.getLog10PosteriorOfAFGT0(), o2.getLog10PosteriorOfAFGT0());
        }
    }

    private final static CompareAFCalcResultsByPNonRef compareAFCalcResultsByPNonRef = new CompareAFCalcResultsByPNonRef();

    final ReferenceDiploidExactAFCalc refModel;

    protected IndependentAllelesDiploidExactAFCalc(int nSamples, int maxAltAlleles, int maxAltAllelesForIndels, final int ploidy) {
        super(nSamples, maxAltAlleles, maxAltAllelesForIndels, ploidy);
        refModel = new ReferenceDiploidExactAFCalc(nSamples, 1, 1, ploidy);
    }

    @Override
    protected StateTracker makeMaxLikelihood(VariantContext vc, AFCalcResultTracker resultTracker) {
        return refModel.makeMaxLikelihood(vc, resultTracker);
    }

    private static class MyAFCalcResult extends AFCalcResult {
        final List<AFCalcResult> supporting;

        private MyAFCalcResult(int[] alleleCountsOfMLE, int nEvaluations, List<Allele> allelesUsedInGenotyping, double[] log10LikelihoodsOfAC, double[] log10PriorsOfAC, Map<Allele, Double> log10pNonRefByAllele, List<AFCalcResult> supporting) {
            super(alleleCountsOfMLE, nEvaluations, allelesUsedInGenotyping, log10LikelihoodsOfAC, log10PriorsOfAC, log10pNonRefByAllele);
            this.supporting = supporting;
        }
    }

    @Override
    public AFCalcResult computeLog10PNonRef(final VariantContext vc,
                                            final double[] log10AlleleFrequencyPriors) {
        final List<AFCalcResult> independentResultTrackers = computeAlleleConditionalExact(vc, log10AlleleFrequencyPriors);
        final List<AFCalcResult> withMultiAllelicPriors = applyMultiAllelicPriors(independentResultTrackers);
        return combineIndependentPNonRefs(vc, withMultiAllelicPriors);
    }


    /**
     *
     * @param vc
     * @param log10AlleleFrequencyPriors
     * @return
     */
    protected List<AFCalcResult> computeAlleleConditionalExact(final VariantContext vc,
                                                               final double[] log10AlleleFrequencyPriors) {
        final List<AFCalcResult> results = new LinkedList<AFCalcResult>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc) ) {
            final AFCalcResult resultTracker = refModel.getLog10PNonRef(subvc, log10AlleleFrequencyPriors);
            results.add(resultTracker);
        }

        return results;
    }

    protected List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;
        final List<VariantContext> vcs = new LinkedList<VariantContext>();

        final List<Allele> afZeroAlleles = new LinkedList<Allele>();
        for ( int altI = 0; altI < nAltAlleles; altI++ ) {
            final Allele altAllele = vc.getAlternateAllele(altI);
            final List<Allele> biallelic = Arrays.asList(vc.getReference(), altAllele);
            vcs.add(biallelicCombinedGLs(vc, biallelic, afZeroAlleles, altI + 1));
            //afZeroAlleles.add(altAllele);
        }

        return vcs;
    }

    protected VariantContext biallelicCombinedGLs(final VariantContext rootVC, final List<Allele> biallelic, final List<Allele> afZeroAlleles, final int allele2) {
        if ( rootVC.isBiallelic() ) {
            if ( ! afZeroAlleles.isEmpty() ) throw new IllegalArgumentException("Root VariantContext is biallelic but afZeroAlleles wasn't empty: " + afZeroAlleles);
            return rootVC;
        } else {
            final Set<Integer> allelesToDiscard = new HashSet<Integer>(rootVC.getAlleleIndices(afZeroAlleles));
            final int nAlts = rootVC.getNAlleles() - 1;
            final List<Genotype> biallelicGenotypes = new ArrayList<Genotype>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLs(g, allele2, allelesToDiscard, nAlts));

            final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
            vcb.alleles(biallelic);
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
     * Supports the additional mode of simply dropping GLs whose allele index occurs in allelesToDiscard.  This is
     * useful in the case where you want to drop alleles (not combine them), such as above:
     *
     * AA AB BB AC BC CC
     *
     * and we want to get the bi-allelic GLs for X/B, where X is everything not B, but dropping C (index 2)
     *
     * XX = AA (since X = A and C is dropped)
     * XB = AB
     * BB = BB
     *
     * This allows us to recover partial GLs the correspond to any allele in allelesToDiscard having strictly
     * AF == 0.
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @Requires({"original.hasLikelihoods()", "! allelesToDiscard.contains(altIndex)"})
    @Ensures({"result.hasLikelihoods()", "result.getPL().length == 3"})
    protected Genotype combineGLs(final Genotype original, final int altIndex, final Set<Integer> allelesToDiscard, final int nAlts ) {
        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(new int[]{0,0,0}).alleles(BIALLELIC_NOCALL).make();

        if ( altIndex < 1 || altIndex > nAlts ) throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);

        final double[] normalizedPr = MathUtils.normalizeFromLog10(GenotypeLikelihoods.fromPLs(original.getPL()).getAsVector());
        final double[] biAllelicPr = new double[3];

        for ( int index = 0; index < normalizedPr.length; index++ ) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(index);

            // just continue if we shouldn't include the pair because it's in the discard set
            if ( discardAllelePair(pair, allelesToDiscard) )
                continue;

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

    protected boolean discardAllelePair(final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair, Set<Integer> allelesToDiscard) {
        return allelesToDiscard.contains(pair.alleleIndex1) || allelesToDiscard.contains(pair.alleleIndex2);
    }

    protected List<AFCalcResult> applyMultiAllelicPriors(final List<AFCalcResult> conditionalPNonRefResults) {
        final ArrayList<AFCalcResult> sorted = new ArrayList<AFCalcResult>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, compareAFCalcResultsByPNonRef);

        final double log10SingleAllelePriorOfAFGt0 = conditionalPNonRefResults.get(0).getLog10PriorOfAFGT0();

        for ( int i = 0; i < sorted.size(); i++ ) {
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
     * @param sortedResultsWithThetaNPriors the pNonRef result for each allele independently
     */
    protected AFCalcResult combineIndependentPNonRefs(final VariantContext vc,
                                                      final List<AFCalcResult> sortedResultsWithThetaNPriors) {
        int nEvaluations = 0;
        final int nAltAlleles = sortedResultsWithThetaNPriors.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final double[] log10PriorsOfAC = new double[2];
        final Map<Allele, Double> log10pNonRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

        // this value is a sum in log space
        double log10PosteriorOfACEq0Sum = 0.0;

        for ( final AFCalcResult sortedResultWithThetaNPriors : sortedResultsWithThetaNPriors ) {
            final Allele altAllele = sortedResultWithThetaNPriors.getAllelesUsedInGenotyping().get(1);
            final int altI = vc.getAlleles().indexOf(altAllele) - 1;

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = sortedResultWithThetaNPriors.getAlleleCountAtMLE(altAllele);

            log10PriorsOfAC[0] += sortedResultWithThetaNPriors.getLog10PriorOfAFEq0();
            log10PriorsOfAC[1] += sortedResultWithThetaNPriors.getLog10PriorOfAFGT0();

            // the AF > 0 case requires us to store the normalized likelihood for later summation
            if ( sortedResultWithThetaNPriors.getLog10PosteriorOfAFGT0() > MIN_LOG10_CONFIDENCE_TO_INCLUDE_ALLELE_IN_POSTERIOR )
                log10PosteriorOfACEq0Sum += sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0();

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            log10pNonRefByAllele.put(altAllele, sortedResultWithThetaNPriors.getLog10PosteriorOfAFGT0());

            // trivial -- update the number of evaluations
            nEvaluations += sortedResultWithThetaNPriors.nEvaluations;
        }

        // In principle, if B_p = x and C_p = y are the probabilities of being poly for alleles B and C,
        // the probability of being poly is (1 - B_p) * (1 - C_p) = (1 - x) * (1 - y).  We want to estimate confidently
        // log10((1 - x) * (1 - y)) which is log10(1 - x) + log10(1 - y).  This sum is log10PosteriorOfACEq0
        final double log10PosteriorOfACGt0 = Math.max(Math.log10(1 - Math.pow(10, log10PosteriorOfACEq0Sum)), MathUtils.LOG10_P_OF_ZERO);
        final double[] log10LikelihoodsOfAC = new double[] {
                // L + prior = posterior => L = poster - prior
                log10PosteriorOfACEq0Sum - log10PriorsOfAC[0],
                log10PosteriorOfACGt0 - log10PriorsOfAC[1]
        };

        return new MyAFCalcResult(alleleCountsOfMLE, nEvaluations, vc.getAlleles(),
                MathUtils.normalizeFromLog10(log10LikelihoodsOfAC, true),   // necessary to ensure all values < 0
                MathUtils.normalizeFromLog10(log10PriorsOfAC, true),        // priors incorporate multiple alt alleles, must be normalized
                log10pNonRefByAllele, sortedResultsWithThetaNPriors);
    }
}
