package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;


/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele)
 * Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class MappingQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities")); }

    protected void fillQualsFromPileup(final List<Allele> allAlleles,
                                       final int refLoc,
                                       final ReadBackedPileup pileup,
                                       final PerReadAlleleLikelihoodMap likelihoodMap,
                                       final List<Double> refQuals, final List<Double> altQuals) {

        if (pileup != null && likelihoodMap == null) {
            // old UG snp-only path through the annotations
            for ( final PileupElement p : pileup ) {
                if ( isUsableBase(p) ) {
                    if ( allAlleles.get(0).equals(Allele.create(p.getBase(), true)) ) {
                        refQuals.add((double)p.getMappingQual());
                    } else if ( allAlleles.contains(Allele.create(p.getBase()))) {
                        altQuals.add((double)p.getMappingQual());
                    }
                }
            }
            return;
        }
        for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : likelihoodMap.getLikelihoodReadMap().entrySet()) {
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            // BUGBUG: There needs to be a comparable isUsableBase check here
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add((double)el.getKey().getMappingQuality());
            else if (allAlleles.contains(a))
                altQuals.add((double)el.getKey().getMappingQuality());
        }
    }

 }