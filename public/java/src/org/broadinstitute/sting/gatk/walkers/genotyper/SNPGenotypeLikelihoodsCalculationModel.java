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
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class SNPGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    // the alternate allele with the largest sum of quality scores
    protected Byte bestAlternateAllele = null;

    private final boolean useAlleleFromVCF;

    protected SNPGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
    }

    public Allele getLikelihoods(RefMetaDataTracker tracker,
                                 ReferenceContext ref,
                                 Map<String, AlignmentContext> contexts,
                                 AlignmentContextUtils.ReadOrientation contextType,
                                 GenotypePriors priors,
                                 Map<String, MultiallelicGenotypeLikelihoods> GLs,
                                 Allele alternateAlleleToUse,
                                 boolean useBAQedPileup) {

        if ( !(priors instanceof DiploidSNPGenotypePriors) )
            throw new StingException("Only diploid-based SNP priors are supported in the SNP GL model");

        byte refBase = ref.getBase();
        Allele refAllele = Allele.create(refBase, true);

        // find the alternate allele with the largest sum of quality scores
        if ( alternateAlleleToUse != null ) {
            bestAlternateAllele = alternateAlleleToUse.getBases()[0];
        } else if ( useAlleleFromVCF ) {
            VariantContext vc = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), true, logger, UAC.alleles);

            // ignore places where we don't have a variant
            if ( vc == null )
                return null;

            if ( !vc.isBiallelic() ) {
                // for multi-allelic sites go back to the reads and find the most likely alternate allele
                initializeBestAlternateAllele(refBase, contexts, useBAQedPileup);
            } else {
                bestAlternateAllele = vc.getAlternateAllele(0).getBases()[0];
            }
        } else {
            initializeBestAlternateAllele(refBase, contexts, useBAQedPileup);
        }

        // if there are no non-ref bases...
        if ( bestAlternateAllele == null ) {
            // if we only want variants, then we don't need to calculate genotype likelihoods
            if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                return refAllele;

            // otherwise, choose any alternate allele (it doesn't really matter)
            bestAlternateAllele = (byte)(refBase != 'A' ? 'A' : 'C');
        }

        Allele altAllele = Allele.create(bestAlternateAllele, false);

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            if( useBAQedPileup ) { pileup = createBAQedPileup( pileup ); }

            // create the GenotypeLikelihoods object
            DiploidSNPGenotypeLikelihoods GL = new DiploidSNPGenotypeLikelihoods((DiploidSNPGenotypePriors)priors, UAC.PCR_error);
            int nGoodBases = GL.add(pileup, true, true, UAC.MIN_BASE_QUALTY_SCORE);
            if ( nGoodBases == 0 )
                continue;

            double[] likelihoods = GL.getLikelihoods();

            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(refBase);
            DiploidGenotype hetGenotype = DiploidGenotype.createDiploidGenotype(refBase, bestAlternateAllele);
            DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(bestAlternateAllele);
            ArrayList<Allele> aList = new ArrayList<Allele>();
            aList.add(refAllele);
            aList.add(altAllele);
            double[] dlike = new double[]{likelihoods[refGenotype.ordinal()],likelihoods[hetGenotype.ordinal()],likelihoods[homGenotype.ordinal()]} ;

            // normalize in log space so that max element is zero.
            GLs.put(sample.getKey(), new MultiallelicGenotypeLikelihoods(sample.getKey(),
                    aList,  MathUtils.normalizeFromLog10(dlike, false, true), getFilteredDepth(pileup)));
        }

        return refAllele;
    }

    protected void initializeBestAlternateAllele(byte ref, Map<String, AlignmentContext> contexts, boolean useBAQedPileup) {
        int[] qualCounts = new int[4];

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // calculate the sum of quality scores for each base
            ReadBackedPileup pileup = useBAQedPileup ? createBAQedPileup( sample.getValue().getBasePileup() ) : sample.getValue().getBasePileup();
            for ( PileupElement p : pileup ) {
                // ignore deletions
                if ( p.isDeletion() || (! p.isReducedRead() && p.getQual() < UAC.MIN_BASE_QUALTY_SCORE ))
                    continue;

                final int index = BaseUtils.simpleBaseToBaseIndex(p.getBase());
                if ( index >= 0 ) {
                    qualCounts[index] += p.getQual();
                }
            }
        }

        // set the non-ref base with maximum quality score sum
        int maxCount = 0;
        bestAlternateAllele = null;
        for ( byte altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;
            int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
            if ( qualCounts[index] > maxCount ) {
                maxCount = qualCounts[index];
                bestAlternateAllele = altAllele;
            }
        }
    }

    public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
    }

    public class BAQedPileupElement extends PileupElement {
        public BAQedPileupElement( final PileupElement PE ) {
            super(PE.getRead(), PE.getOffset());
        }

        @Override
        public byte getQual( final int offset ) { return BAQ.calcBAQFromTag(getRead(), offset, true); }
    }

}