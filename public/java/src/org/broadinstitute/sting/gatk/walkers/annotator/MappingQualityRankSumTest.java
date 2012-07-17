package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;


/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele)
 * Note that the mapping quality rank sum test can not be calculated for homozygous sites.
 */
public class MappingQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities")); }

    protected void fillQualsFromPileup(byte ref, List<Byte> alts, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        for ( final PileupElement p : pileup ) {
            if ( isUsableBase(p) ) {
                if ( p.getBase() == ref ) {
                    refQuals.add((double)p.getMappingQual());
                } else if ( alts.contains(p.getBase()) ) {
                    altQuals.add((double)p.getMappingQual());
                }
            }
        }
    }

    protected void fillQualsFromPileup(final Allele ref, final List<Allele> alts, final int refLoc, final Map<Allele, List<GATKSAMRecord>> stratifiedContext, final List<Double> refQuals, final List<Double> altQuals) {
        for ( final Map.Entry<Allele, List<GATKSAMRecord>> alleleBin : stratifiedContext.entrySet() ) {
            final boolean matchesRef = ref.equals(alleleBin.getKey());
            final boolean matchesAlt = alts.contains(alleleBin.getKey());
            if ( !matchesRef && !matchesAlt )
                continue;

            for ( final GATKSAMRecord read : alleleBin.getValue() ) {
                if ( matchesRef )
                    refQuals.add((double)read.getMappingQuality());
                else
                    altQuals.add((double)read.getMappingQuality());
            }
        }
    }

    protected void fillIndelQualsFromPileup(ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        // equivalent is whether indel likelihoods for reads corresponding to ref allele are more likely than reads corresponding to alt allele ?
        HashMap<PileupElement,LinkedHashMap<Allele,Double>> indelLikelihoodMap = IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap();
        for (final PileupElement p: pileup) {
            if (indelLikelihoodMap.containsKey(p) && p.getMappingQual() != 0 && p.getMappingQual() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
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
                    refQuals.add((double)p.getMappingQual());
                else if (altLikelihood > refLikelihood + INDEL_LIKELIHOOD_THRESH)
                    altQuals.add((double)p.getMappingQual());
            }
        }
    }
    
}