package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.List;
import java.util.Arrays;


public class MappingQualityRankSumTest /*extends RankSumTest*/ {

    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Het Vs. Ref Read Mapping Qualities")); }

    protected void fillQualsFromPileup(byte ref, char alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( PileupElement p : pileup ) {
            // ignore deletions
            if ( p.isDeletion() )
                continue;

            if ( p.getBase() == ref )
                refQuals.add(p.getMappingQual());
            else if ( (char)p.getBase() == alt )
                altQuals.add(p.getMappingQual());
        }
    }
}