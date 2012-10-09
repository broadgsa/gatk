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
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class IndependentAllelesDiploidExactAFCalc extends DiploidExactAFCalc {
    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
    final ReferenceDiploidExactAFCalc refModel;

    public IndependentAllelesDiploidExactAFCalc(final int nSamples, final int maxAltAlleles) {
        super(nSamples, maxAltAlleles);
        refModel = new ReferenceDiploidExactAFCalc(nSamples, 1);
    }

    public IndependentAllelesDiploidExactAFCalc(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        refModel = new ReferenceDiploidExactAFCalc(nSamples, 1);
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
        final List<AFCalcResult> independentResultTrackers = computeLog10PNonRefForEachAllele(vc, log10AlleleFrequencyPriors);
        return combineIndependentPNonRefs(vc, log10LikelihoodOfRef, independentResultTrackers, log10AlleleFrequencyPriors);
    }

    protected final double computelog10LikelihoodOfRef(final VariantContext vc) {
        // this value just the likelihood of AF == 0 in the special constrained multi-allelic calculation
        final List<double[]> allGLs = getGLs(vc.getGenotypes(), false);
        double log10LikelihoodOfHomRef = 0.0;

        // TODO -- can be easily optimized (currently looks at all GLs via getGLs)
        for ( int i = 0; i < allGLs.size(); i++ ) {
            final double[] GLs = allGLs.get(i);
            log10LikelihoodOfHomRef += GLs[0];
        }

        return log10LikelihoodOfHomRef;

//        // this value just the likelihood of AF == 0 in the special constrained multi-allelic calculation
//        final List<double[]> allGLs = getGLs(vc.getGenotypes(), false);
//        final double[] log10LikelihoodOfHomRefs = new double[allGLs.size()];
//
//        // TODO -- can be easily optimized (currently looks at all GLs via getGLs)
//        for ( int i = 0; i < allGLs.size(); i++ ) {
//            final double[] GLs = allGLs.get(i);
//            log10LikelihoodOfHomRefs[i] = GLs[0];
//        }
//
//        return MathUtils.log10sumLog10(log10LikelihoodOfHomRefs);
    }

    protected List<AFCalcResult> computeLog10PNonRefForEachAllele(final VariantContext vc,
                                                                  final double[] log10AlleleFrequencyPriors) {
        final int nAltAlleles = vc.getNAlleles() - 1;
        final List<AFCalcResult> resultTrackers = new ArrayList<AFCalcResult>(nAltAlleles);

        for ( int altI = 0; altI < nAltAlleles; altI++ ) {
            final List<Allele> biallelic = Arrays.asList(vc.getReference(), vc.getAlternateAllele(altI));
            final VariantContext subvc = biallelicCombinedGLs(vc, biallelic, altI + 1);
            final AFCalcResult resultTracker = refModel.getLog10PNonRef(subvc, log10AlleleFrequencyPriors);
            resultTrackers.add(resultTracker);
        }

        return resultTrackers;
    }

    protected VariantContext biallelicCombinedGLs(final VariantContext rootVC, final List<Allele> biallelic, final int allele2) {
        if ( rootVC.isBiallelic() )
            return rootVC;
        else {
            final int nAlts = rootVC.getNAlleles() - 1;
            final List<Genotype> biallelicGenotypes = new ArrayList<Genotype>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLs(g, allele2, nAlts));

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
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @Requires("original.hasLikelihoods()")
    @Ensures({"result.hasLikelihoods()", "result.getPL().length == 3"})
    protected Genotype combineGLs(final Genotype original, final int altIndex, final int nAlts ) {
        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(new int[]{0,0,0}).alleles(BIALLELIC_NOCALL).make();

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

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Takes each independent result and merges it into the final result object
     *
     * Suppose you have L_af=0_1 = -1 and L_af=0_1 = -2 and L_af=1_1 = -3 and L_af=1_2 = 0.  What does this mean?
     * If says that along dimension 1, the AF is more likely to be ref (-1 vs. -3) while along dimension 2
     * you are more likely to be alt (-2 vs. 0).  The question is how to combine these into a meaningful
     * composite likelihood.  What we are interested in is:
     *
     *     L(AF == 0 for all dimensions) vs. L(AF > 0 for any dimension)
     *
     * So what are these quantities?  The problem is that the likelihoods aren't normalized, so we really cannot
     * just add them together.  What we really need are normalized probabilities so that we can compute:
     *
     * P(AF == 0 for all dimensions) => product_i for P(AF == 0, i)
     * P(AF > 0 for any dimension)   => sum_i for P(AF > 0, i)
     *
     * These probabilities can be computed straight off the likelihoods without a prior.  It's just the prior-free
     * normalization of the two likelihoods.
     *
     * @param independentPNonRefs the pNonRef result for each allele independently
     */
    protected AFCalcResult combineIndependentPNonRefs(final VariantContext vc,
                                                      final double log10LikelihoodsOfACEq0,
                                                      final List<AFCalcResult> independentPNonRefs,
                                                      final double[] log10AlleleFrequencyPriors) {
        int nEvaluations = 0;
        final int nAltAlleles = independentPNonRefs.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final double[] log10PriorsOfAC = new double[2];
        final Map<Allele, Double> log10pNonRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

        // this value is a sum in real space so we need to store values to sum up later
        final double[] log10LikelihoodsOfACGt0 = new double[nAltAlleles];

        // TODO -- need to apply theta^alt prior after sorting by MLE

        int altI = 0;
        for ( final AFCalcResult independentPNonRef : independentPNonRefs ) {
            final Allele altAllele = vc.getAlternateAllele(altI);

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = independentPNonRef.getAlleleCountAtMLE(altAllele);

            // TODO -- figure out real value, this is a temp (but good) approximation
            if ( altI == 0 ) {
                log10PriorsOfAC[0] = independentPNonRef.getLog10PriorOfAFEq0();
                log10PriorsOfAC[1] = independentPNonRef.getLog10PriorOfAFGT0();
            }

            // now we effectively have flat prior'd posteriors
            final double[] log10NormalizedLikelihoods = MathUtils.normalizeFromLog10(
                    new double[]{
                            independentPNonRef.getLog10LikelihoodOfAFEq0(),
                            independentPNonRef.getLog10LikelihoodOfAFGT0() },
                    true);

            // the AF > 0 case requires us to store the normalized likelihood for later summation
            log10LikelihoodsOfACGt0[altI] = log10NormalizedLikelihoods[1];

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
                log10PriorsOfAC, log10pNonRefByAllele, independentPNonRefs);
    }
}
