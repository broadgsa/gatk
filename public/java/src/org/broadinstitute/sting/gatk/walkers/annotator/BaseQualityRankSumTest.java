package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;


/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for base qualities (ref bases vs. bases of the alternate allele).
 * Note that the base quality rank sum test can not be calculated for homozygous sites.
 */
public class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation {
    public List<String> getKeyNames() { return Arrays.asList("BaseQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("BaseQRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities")); }

    protected void fillQualsFromPileup(byte ref, List<Byte> alts, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) ) {
                if ( p.getBase() == ref )
                    refQuals.add((double)p.getQual());
                else if ( alts.contains(p.getBase()) )
                    altQuals.add((double)p.getQual());
            }
        }
    }
    protected void fillQualsFromPileup(final Allele ref, final List<Allele> alts, final int refLoc, final Map<Allele, List<GATKSAMRecord>> stratifiedContext, final List<Double> refQuals, final List<Double> altQuals) {
        // TODO -- implement me; how do we pull out the correct offset from the read?
        return;

/*
        for ( final Map.Entry<Allele, List<GATKSAMRecord>> alleleBin : stratifiedContext.entrySet() ) {
            final boolean matchesRef = ref.equals(alleleBin.getKey());
            final boolean matchesAlt = alts.contains(alleleBin.getKey());
            if ( !matchesRef && !matchesAlt )
                continue;

            for ( final GATKSAMRecord read : alleleBin.getValue() ) {

                if ( isUsableBase(p) ) {
                    if ( matchesRef )
                        refQuals.add((double)p.getQual());
                    else
                        altQuals.add((double)p.getQual());
                }
            }
        }
*/
    }

    protected void fillIndelQualsFromPileup(ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        // equivalent is whether indel likelihoods for reads corresponding to ref allele are more likely than reads corresponding to alt allele ?
        HashMap<PileupElement,LinkedHashMap<Allele,Double>> indelLikelihoodMap = IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap();
        for (final PileupElement p: pileup) {
            if (indelLikelihoodMap.containsKey(p)) {
                // retrieve likelihood information corresponding to this read
                LinkedHashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                // by design, first element in LinkedHashMap was ref allele
                double refLikelihood=0.0, altLikelihood=Double.NEGATIVE_INFINITY;

                for (Allele a : el.keySet()) {

                    if (a.isReference())
                        refLikelihood =el.get(a);
                    else {
                        double like = el.get(a);
                        if (like >= altLikelihood)
                            altLikelihood = like;
                    }
                }
                if (refLikelihood > altLikelihood + INDEL_LIKELIHOOD_THRESH)
                    refQuals.add(-10.0*refLikelihood);
                else if (altLikelihood > refLikelihood + INDEL_LIKELIHOOD_THRESH)
                    altQuals.add(-10.0*altLikelihood);
            }
        }
    }

}