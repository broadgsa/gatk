package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/30/11
 */

public class ReadPosRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("ReadPosRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Alt Vs. Ref read position bias")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) ) {
                int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p.getOffset(), 0, 0);
                final int numAlignedBases = AlignmentUtils.getNumAlignedBases(p.getRead());
                if( readPos > numAlignedBases / 2 ) {
                    readPos = numAlignedBases - ( readPos + 1 );
                }

                if ( p.getBase() == ref ) {
                    refQuals.add( readPos );
                } else if ( p.getBase() == alt ) {
                    altQuals.add( readPos );
                }
            }
        }
    }
}
