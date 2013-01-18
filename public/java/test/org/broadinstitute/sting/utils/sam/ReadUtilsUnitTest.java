/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
        int myStart, boundary;

        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, cigar);
        read.setMateAlignmentStart(mateStart);
        read.setInferredInsertSize(fragmentSize);

        // Test case 1: positive strand, first read
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, myStart + fragmentSize + 1);

        // Test case 2: positive strand, second read
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, myStart + fragmentSize + 1);

        // Test case 3: negative strand, second read
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 4: negative strand, first read
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 5: mate is mapped to another chromosome (test both strands)
        read.setInferredInsertSize(0);
        read.setReadNegativeStrandFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setInferredInsertSize(10);

        // Test case 6: read is unmapped
        read.setReadUnmappedFlag(true);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
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
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // second read:
        myStart = 1000;
        read.setAlignmentStart(myStart);
        read.setMateAlignmentStart(980);
        read.setReadNegativeStrandFlag(false);
        boundary = ReadUtils.getAdaptorBoundary(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
    }
}
