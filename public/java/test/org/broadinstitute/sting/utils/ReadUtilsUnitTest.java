package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;


public class ReadUtilsUnitTest extends BaseTest {
    SAMRecord read, reducedRead;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";
    final private static int REDUCED_READ_QUAL = 20;

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
        reducedRead.setAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG, REDUCED_READ_QUAL);
    }

    private void testReadBasesAndQuals(SAMRecord read, int expectedStart, int expectedStop) {
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
        Assert.assertFalse(ReadUtils.isReducedRead(read), "isReducedRead is false for normal read");
        Assert.assertEquals(ReadUtils.getReducedReadQualityTagValue(read), null, "No reduced read tag in normal read");

        Assert.assertTrue(ReadUtils.isReducedRead(reducedRead), "isReducedRead is true for reduced read");
        Assert.assertEquals((int) ReadUtils.getReducedReadQualityTagValue(reducedRead), REDUCED_READ_QUAL, "Reduced read tag is set to expected value");
    }

    @Test
    public void testreducedReadWithReducedQualsWithReducedRead() {
        SAMRecord replacedRead = ReadUtils.reducedReadWithReducedQuals(reducedRead);
        Assert.assertEquals(replacedRead.getReadBases(), reducedRead.getReadBases());
        Assert.assertEquals(replacedRead.getBaseQualities().length, reducedRead.getBaseQualities().length);
        for ( int i = 0; i < replacedRead.getBaseQualities().length; i++)
            Assert.assertEquals(replacedRead.getBaseQualities()[i], REDUCED_READ_QUAL);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testreducedReadWithReducedQualsWithNormalRead() {
        ReadUtils.reducedReadWithReducedQuals(read);
    }

    @Test
    public void testReducedReadPileupElement() {
        PileupElement readp = new PileupElement(read,0);
        PileupElement reducedreadp = new PileupElement(reducedRead,0);

        Assert.assertFalse(readp.isReducedRead());

        Assert.assertTrue(reducedreadp.isReducedRead());
        Assert.assertEquals(reducedreadp.getReducedCount(), 0);
        Assert.assertEquals(reducedreadp.getReducedQual(), REDUCED_READ_QUAL);

    }
}
