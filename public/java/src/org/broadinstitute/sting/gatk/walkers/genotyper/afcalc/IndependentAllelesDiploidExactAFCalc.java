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

public class IndependentAllelesDiploidExactAFCalc extends DiploidExactAFCalc {
    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
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
        final double log10LikelihoodOfRef = computelog10LikelihoodOfRef(vc);
        final List<AFCalcResult> independentResultTrackers = computeAlleleConditionalExact(vc, log10AlleleFrequencyPriors);
        return combineIndependentPNonRefs(vc, log10LikelihoodOfRef, independentResultTrackers, log10AlleleFrequencyPriors);
    }

    protected final double computelog10LikelihoodOfRef(final VariantContext vc) {
        // this value just the likelihood of AF == 0 in the special constrained multi-allelic calculation
        final List<double[]> allGLs = getGLs(vc.getGenotypes(), false);
        double log10LikelihoodOfHomRef = 0.0;

        // TODO -- can be easily optimized (currently looks at all GLs via getGLs)
        for ( int i = 0; i < allGLs.size(); i++ ) {
            final double[] GLs = allGLs.get(i);
            log10LikelihoodOfHomRef += MathUtils.normalizeFromLog10(GLs, true)[0];
        }

        return log10LikelihoodOfHomRef;
    }

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
     * The quantity (1) is the same of all cells except those with AF_c == 0, while (2) is the
     * band at the top where AF_t > 0 and AF_c == 0
     *
     * So (4) is actually (1) + (2).
     *
     * (3) is the direct inverse of the (1) and (2), as we are simultaneously calculating
     *
     * (1*) P(D | AF_c == 0 && AF_t == *) [i.e., T can be anything]
     * (2*) P(D | AF_t == 0 && AF_c == 0) [TODO -- note this value looks like the thing we are supposed to use]
     *
     * This function implements the conditional likelihoods summation for any number of alt
     * alleles (not just the tri-allelic case), where each subsequent variant context is
     * further constrained such that each already considered allele x has AF_x == 0 in the
     * compute.
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
            afZeroAlleles.add(altAllele);
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

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * @param conditionalPNonRefResults the pNonRef result for each allele independently
     */
    protected AFCalcResult combineIndependentPNonRefs(final VariantContext vc,
                                                      final double log10LikelihoodsOfACEq0,
                                                      final List<AFCalcResult> conditionalPNonRefResults,
                                                      final double[] log10AlleleFrequencyPriors) {
        int nEvaluations = 0;
        final int nAltAlleles = conditionalPNonRefResults.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final double[] log10PriorsOfAC = new double[2];
        final Map<Allele, Double> log10pNonRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

        // this value is a sum in real space so we need to store values to sum up later
        final double[] log10LikelihoodsOfACGt0 = new double[nAltAlleles];
        //double log10LikelihoodsOfACEq0 = 0.0;

        // TODO -- need to apply theta^alt prior after sorting by MLE

        int altI = 0;
        for ( final AFCalcResult independentPNonRef : conditionalPNonRefResults ) {
            final Allele altAllele = vc.getAlternateAllele(altI);

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = independentPNonRef.getAlleleCountAtMLE(altAllele);

            // TODO -- figure out real value, this is a temp (but good) approximation
            if ( altI == 0 ) {
                log10PriorsOfAC[0] = independentPNonRef.getLog10PriorOfAFEq0();
                log10PriorsOfAC[1] = independentPNonRef.getLog10PriorOfAFGT0();
            }

            // the AF > 0 case requires us to store the normalized likelihood for later summation
            //log10LikelihoodsOfACEq0 += independentPNonRef.getLog10LikelihoodOfAFEq0();
            log10LikelihoodsOfACGt0[altI] = independentPNonRef.getLog10LikelihoodOfAFGT0();

            // bind pNonRef for allele to the posterior value of the AF > 0
            // TODO -- should incorporate the theta^alt prior here from the likelihood itself
            log10pNonRefByAllele.put(altAllele, independentPNonRef.getLog10PosteriorOfAFGt0ForAllele(altAllele));

            // trivial -- update the number of evaluations
            nEvaluations += independentPNonRef.nEvaluations;
            altI++;
        }

        // the log10 likelihoods are the sum of the log10 likelihoods across all alt alleles
        final double[] log10LikelihoodsOfAC = new double[]{
                log10LikelihoodsOfACEq0,
                MathUtils.log10sumLog10(log10LikelihoodsOfACGt0)};

        return new MyAFCalcResult(alleleCountsOfMLE, nEvaluations, vc.getAlleles(),
                MathUtils.normalizeFromLog10(log10LikelihoodsOfAC, true, true), // necessary to ensure all values < 0
                log10PriorsOfAC, log10pNonRefByAllele, conditionalPNonRefResults);
    }
}
