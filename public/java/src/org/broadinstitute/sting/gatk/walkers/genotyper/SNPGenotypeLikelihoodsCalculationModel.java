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
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class SNPGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    private static final int MIN_QUAL_SUM_FOR_ALT_ALLELE = 50;

    private boolean ALLOW_MULTIPLE_ALLELES;

    private final boolean useAlleleFromVCF;

    protected SNPGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        ALLOW_MULTIPLE_ALLELES = UAC.MULTI_ALLELIC;
        useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
    }

    public VariantContext getLikelihoods(RefMetaDataTracker tracker,
                                         ReferenceContext ref,
                                         Map<String, AlignmentContext> contexts,
                                         AlignmentContextUtils.ReadOrientation contextType,
                                         GenotypePriors priors,
                                         Allele alternateAlleleToUse,
                                         boolean useBAQedPileup) {

        if ( !(priors instanceof DiploidSNPGenotypePriors) )
            throw new StingException("Only diploid-based SNP priors are supported in the SNP GL model");

        final boolean[] basesToUse = new boolean[4];
        final byte refBase = ref.getBase();
        final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(refBase);

        // start making the VariantContext
        final GenomeLoc loc = ref.getLocus();
        final List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(refBase, true));
        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), alleles);

        // find the alternate allele(s) that we should be using
        if ( alternateAlleleToUse != null ) {
            basesToUse[BaseUtils.simpleBaseToBaseIndex(alternateAlleleToUse.getBases()[0])] = true;
        } else if ( useAlleleFromVCF ) {
            final VariantContext vc = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), true, logger, UAC.alleles);

            // ignore places where we don't have a SNP
            if ( vc == null || !vc.isSNP() )
                return null;

            for ( Allele allele : vc.getAlternateAlleles() )
                basesToUse[BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0])] = true;
        } else {

            determineAlternateAlleles(basesToUse, refBase, contexts, useBAQedPileup);

            // how many alternate alleles are we using?
            int alleleCounter = Utils.countSetBits(basesToUse);

            // if there are no non-ref alleles...
            if ( alleleCounter == 0 ) {
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    return builder.make();

                // otherwise, choose any alternate allele (it doesn't really matter)
                basesToUse[indexOfRefBase == 0 ? 1 : 0] = true;
             }
        }

        // create the alternate alleles and the allele ordering (the ordering is crucial for the GLs)
        final int numAltAlleles = Utils.countSetBits(basesToUse);
        final int[] alleleOrdering = new int[numAltAlleles + 1];
        alleleOrdering[0] = indexOfRefBase;
        int alleleOrderingIndex = 1;
        int numLikelihoods = 1;
        for ( int i = 0; i < 4; i++ ) {
            if ( i != indexOfRefBase && basesToUse[i] ) {
                alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(i), false));
                alleleOrdering[alleleOrderingIndex++] = i;
                numLikelihoods += alleleOrderingIndex;
            }
        }
        builder.alleles(alleles);

        // create the genotypes; no-call everyone for now
        GenotypesContext genotypes = GenotypesContext.create();
        final List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            if ( useBAQedPileup )
                pileup = createBAQedPileup( pileup );

            // create the GenotypeLikelihoods object
            final DiploidSNPGenotypeLikelihoods GL = new DiploidSNPGenotypeLikelihoods((DiploidSNPGenotypePriors)priors, UAC.PCR_error);
            final int nGoodBases = GL.add(pileup, true, true, UAC.MIN_BASE_QUALTY_SCORE);
            if ( nGoodBases == 0 )
                continue;

            final double[] allLikelihoods = GL.getLikelihoods();
            final double[] myLikelihoods = new double[numLikelihoods];

            int myLikelihoodsIndex = 0;
            for ( int i = 0; i <= numAltAlleles; i++ ) {
                for ( int j = i; j <= numAltAlleles; j++ ) {
                    myLikelihoods[myLikelihoodsIndex++] = allLikelihoods[DiploidGenotype.createDiploidGenotype(alleleOrdering[i], alleleOrdering[j]).ordinal()];
                }
            }

            // normalize in log space so that max element is zero.
            GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(myLikelihoods, false, true));

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.DEPTH_KEY, getFilteredDepth(pileup));
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);
            genotypes.add(new Genotype(sample.getKey(), noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }

        return builder.genotypes(genotypes).make();
    }

    // fills in the allelesToUse array
    protected void determineAlternateAlleles(boolean[] allelesToUse, byte ref, Map<String, AlignmentContext> contexts, boolean useBAQedPileup) {
        int[] qualCounts = new int[4];

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // calculate the sum of quality scores for each base
            ReadBackedPileup pileup = useBAQedPileup ? createBAQedPileup( sample.getValue().getBasePileup() ) : sample.getValue().getBasePileup();
            for ( PileupElement p : pileup ) {
                // ignore deletions
                if ( p.isDeletion() || (!p.isReducedRead() && p.getQual() < UAC.MIN_BASE_QUALTY_SCORE) )
                    continue;

                final int index = BaseUtils.simpleBaseToBaseIndex(p.getBase());
                if ( index >= 0 ) {
                    qualCounts[index] += p.getQual();
                }
            }
        }

        if ( ALLOW_MULTIPLE_ALLELES ) {
            for ( byte altAllele : BaseUtils.BASES ) {
                if ( altAllele == ref )
                    continue;
                int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
                if ( qualCounts[index] >= MIN_QUAL_SUM_FOR_ALT_ALLELE ) {
                    allelesToUse[index] = true;
                }
            }
        } else {
            // set the non-ref base which has the maximum quality score sum
            int maxCount = 0;
            int indexOfMax = 0;
            for ( byte altAllele : BaseUtils.BASES ) {
                if ( altAllele == ref )
                    continue;
                int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
                if ( qualCounts[index] > maxCount ) {
                    maxCount = qualCounts[index];
                    indexOfMax = index;
                }
            }

            if ( maxCount > 0 )
                allelesToUse[indexOfMax] = true;
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