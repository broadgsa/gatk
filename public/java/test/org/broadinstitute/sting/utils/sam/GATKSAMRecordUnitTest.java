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

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


public class GATKSAMRecordUnitTest extends BaseTest {
    GATKSAMRecord read, reducedRead;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";
    final private static byte[] REDUCED_READ_COUNTS = new byte[]{10, 20, 30, 40, 1};
    final private static byte[] REDUCED_READ_COUNTS_TAG = new byte[]{10, 10, 20, 30, -9};  // just the offsets

    @BeforeClass
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
        PileupElement readp = LocusIteratorByState.createPileupForReadAndOffset(read, 0);
        PileupElement reducedreadp = LocusIteratorByState.createPileupForReadAndOffset(reducedRead, 0);

        Assert.assertFalse(readp.getRead().isReducedRead());

        Assert.assertTrue(reducedreadp.getRead().isReducedRead());
        Assert.assertEquals(reducedreadp.getRepresentativeCount(), REDUCED_READ_COUNTS[0]);
        Assert.assertEquals(reducedreadp.getQual(), readp.getQual());
    }

    @Test
    public void testGetOriginalAlignments() {
        final byte [] bases = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals = {20 , 20 , 20 , 20 , 20 , 20 , 20 , 20 };
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, "6M");

        // A regular read with all matches
        Assert.assertEquals(read.getAlignmentStart(), read.getOriginalAlignmentStart());
        Assert.assertEquals(read.getAlignmentEnd(), read.getOriginalAlignmentEnd());

        // Alignment start shifted
        int alignmentShift = 2;
        read.setAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT, alignmentShift);
        Assert.assertEquals(read.getAlignmentStart() + alignmentShift, read.getOriginalAlignmentStart());
        Assert.assertEquals(read.getAlignmentEnd(), read.getOriginalAlignmentEnd());

        // Both alignments shifted
        read.setAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT, alignmentShift);
        Assert.assertEquals(read.getAlignmentStart() + alignmentShift, read.getOriginalAlignmentStart());
        Assert.assertEquals(read.getAlignmentEnd() - alignmentShift, read.getOriginalAlignmentEnd());

        // Alignment end shifted
        read.setAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT, null);
        Assert.assertEquals(read.getAlignmentStart(), read.getOriginalAlignmentStart());
        Assert.assertEquals(read.getAlignmentEnd() - alignmentShift, read.getOriginalAlignmentEnd());
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

    @Test
    public void testGetReducedCountsIsCorrect() {
        final byte[] counts = reducedRead.getReducedReadCounts();
        Assert.assertNotSame(counts, reducedRead.getAttribute(GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG));
        for ( int i = 0; i < counts.length; i++ )
            Assert.assertEquals(counts[i], reducedRead.getReducedCount(i), "Reduced counts vector not equal to getReducedCount(i) at " + i);
    }
}
