package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.List;
import java.util.Arrays;


public class MappingQualityRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("MQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MQRankSum", 1, VCFHeaderLineType.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Alt Vs. Ref read mapping qualities")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) && p.getMappingQual() < 254 ) { // 254 and 255 are special mapping qualities used as a code by aligners
                if ( p.getBase() == ref ) {
                    refQuals.add(p.getMappingQual());
                } else if ( p.getBase() == alt ) {
                    altQuals.add(p.getMappingQual());
                }
            }
        }
    }
}