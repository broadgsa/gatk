package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.List;


public class MappingQualityRankSumTest extends RankSumTest {

    public String getKeyName() { return "MQRankSum"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("MQRankSum", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Het Vs. Ref Read Mapping Qualities"); }

    protected void fillQualsFromPileup(char ref, char alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( PileupElement p : pileup ) {
            // ignore deletions
            if ( p.isDeletion() )
                continue;

            char base = (char)p.getBase();
            if ( base == ref )
                refQuals.add(p.getMappingQual());
            else if ( base == alt )
                altQuals.add(p.getMappingQual());
        }
    }
}