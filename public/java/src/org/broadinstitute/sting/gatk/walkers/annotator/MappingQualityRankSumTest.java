package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;


/**
 * The phred-scaled p-value (u-based z-approximation) from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele)
 */
public class MappingQualityRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        for ( final PileupElement p : pileup ) {
            if ( isUsableBase(p) ) {
                if ( p.getBase() == ref ) {
                    refQuals.add((double)p.getMappingQual());
                } else if ( p.getBase() == alt ) {
                    altQuals.add((double)p.getMappingQual());
                }
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