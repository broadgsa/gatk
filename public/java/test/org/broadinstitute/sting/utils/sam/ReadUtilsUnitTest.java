package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;


public class ReadUtilsUnitTest extends BaseTest {
    GATKSAMRecord read, reducedRead;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";
    final private static byte[] REDUCED_READ_COUNTS = new byte[]{10, 20, 30, 40, 1};
    final private static byte[] REDUCED_READ_COUNTS_TAG = new byte[]{10, 10, 20, 30, -9};  // just the offsets

    @BeforeTest
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));

        reducedRead = ArtificialSAMUtils.createArtificialRead(header, "reducedRead", 0, 1, BASES.length());
        reducedRead.setReadBases(BASES.getBytes());
        reducedRead.setBaseQualityString(QUALS);
        reducedRead.setAttribute(GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG, REDUCED_READ_COUNTS_TAG);
    }

    @Test
    public void testReducedReads() {
        Assert.assertFalse(read.isReducedRead(), "isReducedRead is false for normal read");
        Assert.assertEquals(read.getReducedReadCounts(), null, "No reduced read tag in normal read");

        Assert.assertTrue(reducedRead.isReducedRead(), "isReducedRead is true for reduced read");
        for (int i = 0; i < reducedRead.getReadLength(); i++) {
            Assert.assertEquals(reducedRead.getReducedCount(i), REDUCED_READ_COUNTS[i], "Reduced read count not set to the expected value at " + i);
        }
    }

    @Test
    public void testReducedReadPileupElement() {
        PileupElement readp = new PileupElement(read, 0);
        PileupElement reducedreadp = new PileupElement(reducedRead, 0);

        Assert.assertFalse(readp.isReducedRead());

        Assert.assertTrue(reducedreadp.isReducedRead());
        Assert.assertEquals(reducedreadp.getRepresentativeCount(), REDUCED_READ_COUNTS[0]);
        Assert.assertEquals(reducedreadp.getQual(), readp.getQual());
    }

    @Test
    public void testGetAdaptorBoundary() {
        final byte[] bases = {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
        final byte[] quals = {30, 30, 30, 30, 30, 30, 30, 30};
        final String cigar = "8M";
        final int fragmentSize = 10;
        final int mateStart = 1000;
        final int BEFORE = mateStart - 2;
        final int AFTER = mateStart + 2;
        Integer myStart, boundary;

        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, cigar);
        read.setMateAlignmentStart(mateStart);
        read.setInferredInsertSize(fragmentSize);

        // Test case 1: positive strand, first read
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary.intValue(), myStart + fragmentSize + 1);

        // Test case 2: positive strand, second read
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary.intValue(), myStart + fragmentSize + 1);

        // Test case 3: negative strand, second read
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary.intValue(), mateStart - 1);

        // Test case 4: negative strand, first read
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary.intValue(), mateStart - 1);

        // Test case 5: mate is mapped to another chromosome (test both strands)
        read.setInferredInsertSize(0);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertNull(boundary);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertNull(boundary);
        read.setInferredInsertSize(10);

        // Test case 6: read is unmapped
        read.setReadUnmappedFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertNull(boundary);
        read.setReadUnmappedFlag(false);

        // Test case 7:  reads don't overlap and look like this:
        //    <--------|
        //                 |------>
        // first read:
        myStart = 980;
        read.setAlignmentStart(myStart);
        read.setInferredInsertSize(20);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertNull(boundary);

        // second read:
        myStart = 1000;
        read.setAlignmentStart(myStart);
        read.setMateAlignmentStart(980);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertNull(boundary);
    }
}
