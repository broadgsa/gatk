package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.Cigar;
import org.broadinstitute.sting.utils.clipping.ReadClipperTestUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.List;

/**
 * Checks that the Bad Cigar filter works for all kinds of wonky cigars
 *
 * @author Mauricio Carneiro
 * @since 3/20/12
 */
public class BadCigarFilterUnitTest {

    public static final String[] BAD_CIGAR_LIST = {
            "2D4M",               // starting with multiple deletions
            "4M2D",               // ending with multiple deletions
            "3M1I1D",             // adjacent indels AND ends in deletion
            "1M1I1D2M",           // adjacent indels I->D
            "1M1D2I1M",           // adjacent indels D->I
            "1M1I2M1D",           // ends in single deletion with insertion in the middle
            "4M1D",               // ends in single deletion
            "1D4M",               // starts with single deletion
            "2M1D1D2M",           // adjacent D's
            "1M1I1I1M",           // adjacent I's
            "1H1D4M",             // starting with deletion after H
            "1S1D3M",             // starting with deletion after S
            "1H1S1D3M",           // starting with deletion after HS
            "4M1D1H",             // ending with deletion before H
            "3M1D1S",             // ending with deletion before S
            "3M1D1S1H",           // ending with deletion before HS
            "10M2H10M",           // H in the middle
            "10M2S10M",           // S in the middle
            "1H1S10M2S10M1S1H",    // deceiving S in the middle
            "1H1S10M2H10M1S1H"    // deceiving H in the middle
    };

    BadCigarFilter filter;

    @BeforeClass
    public void init() {
        filter = new BadCigarFilter();
    }

    @Test(enabled = true)
    public void testWonkyCigars () {
        for (String cigarString : BAD_CIGAR_LIST) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
            Assert.assertTrue(filter.filterOut(read), read.getCigarString());
        }
    }

    @Test(enabled = true)
    public void testGoodCigars() {
        List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(10);
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            Assert.assertFalse(filter.filterOut(read), read.getCigarString());
        }
    }
}
