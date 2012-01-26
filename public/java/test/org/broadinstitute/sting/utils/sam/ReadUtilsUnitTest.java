package org.broadinstitute.sting.utils.sam;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;


public class ReadUtilsUnitTest extends BaseTest {
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
