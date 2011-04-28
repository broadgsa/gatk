package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.List;
import java.util.Arrays;


public class BaseQualityRankSumTest extends RankSumTest {
    public List<String> getKeyNames() { return Arrays.asList("BaseQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("BaseQRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) ) {
                if ( p.getBase() == ref ) {
                    refQuals.add((int)p.getQual());
                } else if ( p.getBase() == alt ) {
                    altQuals.add((int)p.getQual());
                }
            }
        }
    }
}