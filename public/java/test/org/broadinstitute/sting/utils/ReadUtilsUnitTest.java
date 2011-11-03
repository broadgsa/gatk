package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;


public class ReadUtilsUnitTest extends BaseTest {
    GATKSAMRecord read, reducedRead;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";
    final private static byte[] REDUCED_READ_COUNTS = new byte[]{10, 20, 30, 40};

    @BeforeTest
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));

        reducedRead = ArtificialSAMUtils.createArtificialRead(header, "reducedRead", 0, 1, BASES.length());
        reducedRead.setReadBases(BASES.getBytes());
        reducedRead.setBaseQualityString(QUALS);
        reducedRead.setAttribute(GATKSAMRecord.REDUCED_READ_QUALITY_TAG, REDUCED_READ_COUNTS);
    }

    private void testReadBasesAndQuals(GATKSAMRecord read, int expectedStart, int expectedStop) {
        SAMRecord clipped = ReadUtils.hardClipBases(read, expectedStart, expectedStop - 1, null);
        String expectedBases = BASES.substring(expectedStart, expectedStop);
        String expectedQuals = QUALS.substring(expectedStart, expectedStop);
        Assert.assertEquals(clipped.getReadBases(), expectedBases.getBytes(), "Clipped bases not those expected");
        Assert.assertEquals(clipped.getBaseQualityString(), expectedQuals, "Clipped quals not those expected");
    }

    @Test public void testNoClip() { testReadBasesAndQuals(read, 0, 4); }
    @Test public void testClip1Front() { testReadBasesAndQuals(read, 1, 4); }
    @Test public void testClip2Front() { testReadBasesAndQuals(read, 2, 4); }
    @Test public void testClip1Back() { testReadBasesAndQuals(read, 0, 3); }
    @Test public void testClip2Back() { testReadBasesAndQuals(read, 0, 2); }

    @Test
    public void testReducedReads() {
        Assert.assertFalse(read.isReducedRead(), "isReducedRead is false for normal read");
        Assert.assertEquals(read.getReducedReadCounts(), null, "No reduced read tag in normal read");

        Assert.assertTrue(reducedRead.isReducedRead(), "isReducedRead is true for reduced read");
        for ( int i = 0; i < reducedRead.getReadLength(); i++) {
            Assert.assertEquals(reducedRead.getReducedCount(i), REDUCED_READ_COUNTS[i], "Reduced read count not set to the expected value at " + i);
        }
    }

    @Test
    public void testReducedReadPileupElement() {
        PileupElement readp = new PileupElement(read,0);
        PileupElement reducedreadp = new PileupElement(reducedRead,0);

        Assert.assertFalse(readp.isReducedRead());

        Assert.assertTrue(reducedreadp.isReducedRead());
        Assert.assertEquals(reducedreadp.getRepresentativeCount(), REDUCED_READ_COUNTS[0]);
        Assert.assertEquals(reducedreadp.getQual(), readp.getQual());
    }
}
