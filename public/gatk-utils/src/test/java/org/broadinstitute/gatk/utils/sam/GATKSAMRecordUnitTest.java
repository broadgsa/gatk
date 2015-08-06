/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


public class GATKSAMRecordUnitTest extends BaseTest {
    GATKSAMRecord read;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";

    @BeforeClass
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));
    }

    @Test
    public void testStrandlessReads() {
        final byte [] bases = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals = {20 , 20 , 20 , 20 , 20 , 20 , 20 , 20 };
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, "6M");
        Assert.assertEquals(read.isStrandless(), false);

        read.setReadNegativeStrandFlag(false);
        Assert.assertEquals(read.isStrandless(), false);
        Assert.assertEquals(read.getReadNegativeStrandFlag(), false);

        read.setReadNegativeStrandFlag(true);
        Assert.assertEquals(read.isStrandless(), false);
        Assert.assertEquals(read.getReadNegativeStrandFlag(), true);

        read.setReadNegativeStrandFlag(true);
        read.setIsStrandless(true);
        Assert.assertEquals(read.isStrandless(), true);
        Assert.assertEquals(read.getReadNegativeStrandFlag(), false, "negative strand flag should return false even through its set for a strandless read");
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testStrandlessReadsFailSetStrand() {
        final byte [] bases = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals = {20 , 20 , 20 , 20 , 20 , 20 , 20 , 20 };
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, "6M");
        read.setIsStrandless(true);
        read.setReadNegativeStrandFlag(true);
    }
}
