package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.List;


public class BaseQualityRankSumTest extends RankSumTest {

    public String getKeyName() { return "BaseQRankSum"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("BaseQRankSum", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Het Vs. Ref Base Qualities"); }

    protected void fillQualsFromPileup(char ref, char alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( PileupElement p : pileup ) {
            // ignore deletions
            if ( p.isDeletion() )
                continue;

            char base = (char)p.getBase();
            if ( base == ref )
                refQuals.add((int)p.getQual());
            else if ( base == alt )
                altQuals.add((int)p.getQual());
        }
    }
}