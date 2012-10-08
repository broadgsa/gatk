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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

    @Override
    public void computeLog10PNonRef(final VariantContext vc,
                                    final double[] log10AlleleFrequencyPriors,
                                    final AFCalcResultTracker resultTracker) {
        refModel.computeLog10PNonRef(vc, log10AlleleFrequencyPriors, resultTracker);
//        final List<AFCalcResult> independentResultTrackers = computeLog10PNonRefForEachAllele(vc, log10AlleleFrequencyPriors);
//        combineIndependentPNonRefs(vc, independentResultTrackers, log10AlleleFrequencyPriors, resultTracker);
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
     * AA AB BB AC BC CC => AA AB+BC CC when altIndex == 1 and nAlts == 2
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @Requires("original.hasLikelihoods()")
    @Ensures({"result.hasLikelihoods()", "result.getPL().length == 3"})
    protected Genotype combineGLs(final Genotype original, final int altIndex, final int nAlts ) {
        if ( altIndex < 1 || altIndex > nAlts ) throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);

        final double[] normalizedPr = MathUtils.normalizeFromLog10(GenotypeLikelihoods.fromPLs(original.getPL()).getAsVector());
        final double[] biAllelicPr = new double[3];
        biAllelicPr[0] = normalizedPr[GenotypeLikelihoods.calculatePLindex(0, 0)];

        for ( int allele1 = 0; allele1 < nAlts+1; allele1++ ) {
            if ( allele1 != altIndex ) {
                final int i = Math.min(altIndex, allele1);
                final int j = Math.max(altIndex, allele1);
                biAllelicPr[1] += normalizedPr[GenotypeLikelihoods.calculatePLindex(i, j)];
            }
        }

        biAllelicPr[2] = normalizedPr[GenotypeLikelihoods.calculatePLindex(altIndex, altIndex)];

        final double[] GLs = new double[3];
        for ( int i = 0; i < GLs.length; i++ ) GLs[i] = Math.log10(biAllelicPr[i]);

        return new GenotypeBuilder(original).PL(GLs).alleles(BIALLELIC_NOCALL).make();
    }

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Takes each independent result and merges it into the final result object
     *
     * @param independentPNonRefs the pNonRef result for each allele independently
     * @param resultTracker the destination for the combined result
     */
    protected void combineIndependentPNonRefs(final VariantContext vc,
                                              final List<AFCalcResult> independentPNonRefs,
                                              final double[] log10AlleleFrequencyPriors,
                                              final AFCalcResultTracker resultTracker) {
//        final int nChrom = vc.getNSamples() * 2;
//
//        resultTracker.reset();
//
//        // both the likelihood and the posterior of AF=0 are the same for all alleles
//        // TODO -- check and ensure this is true
//        resultTracker.setLog10LikelihoodOfAFzero(independentPNonRefs.get(0).getLog10LikelihoodOfAFzero());
//        resultTracker.setLog10PosteriorOfAFzero(independentPNonRefs.get(0).getLog10PosteriorOfAFzero());
//        resultTracker.log10PosteriorMatrixSum = 0.0;
//
//        int altI = 0;
//        for ( final AFCalcResult independentPNonRef : independentPNonRefs ) {
//            resultTracker.log10MLE += independentPNonRef.getLog10MLE();
//
//            // TODO -- technically double counting some posterior mass
//            resultTracker.log10MAP += independentPNonRef.getLog10MAP();
//
//            // TODO -- technically double counting some posterior mass
//            resultTracker.log10PosteriorMatrixSum += independentPNonRef.getLog10PosteriorsMatrixSumWithoutAFzero();
//
//            resultTracker.getAlleleCountsOfMAP()[altI] = independentPNonRef.getAlleleCountsOfMAP()[0];
//            resultTracker.getAlleleCountsOfMLE()[altI] = independentPNonRef.getAlleleCountsOfMLE()[0];
//
//            resultTracker.nEvaluations += independentPNonRef.nEvaluations;
//            altI++;
//        }
    }
}
