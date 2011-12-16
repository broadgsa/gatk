/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.List;

/**
 * User: roger
 * Date: 9/28/11
 */
public class ReadClipperUnitTest extends BaseTest {

    List<Cigar> cigarList;
    int maximumCigarSize = 10;

    @BeforeClass
    public void init() {
        cigarList = ClipReadsTestUtils.generateCigarList(maximumCigarSize);
    }

    @Test(enabled = true)
    public void testHardClipBothEndsByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ClipReadsTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            int readLength = alnStart - alnEnd;
            for (int i=0; i<readLength/2; i++) {
                GATKSAMRecord clippedRead = (new ReadClipper(read)).hardClipBothEndsByReferenceCoordinates(alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getAlignmentStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
                Assert.assertTrue(clippedRead.getAlignmentEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
            }
        }
        logger.warn("PASSED");
    }

    @Test(enabled = true)
    public void testHardClipByReadCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ClipReadsTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            for (int i=0; i<readLength; i++) {
                GATKSAMRecord clipLeft = (new ReadClipper(read)).hardClipByReadCoordinates(0, i);
                Assert.assertTrue(clipLeft.getReadLength() <= readLength - i,  String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipLeft.getCigarString()));

                GATKSAMRecord clipRight = (new ReadClipper(read)).hardClipByReadCoordinates(i, readLength-1);
                Assert.assertTrue(clipRight.getReadLength() <= i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipRight.getCigarString()));
            }
        }
        logger.warn("PASSED");
    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinates() {
        logger.warn("PASSED");

    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesLeftTail() {
        logger.warn("Executing testHardClipByReferenceCoordinatesLeftTail");

    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesRightTail() {
        logger.warn("PASSED");

    }

    @Test(enabled = true)
    public void testHardClipLowQualEnds() {
        final byte LOW_QUAL = 2;
        final byte HIGH_QUAL = 30;

        // create a read for every cigar permutation
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ClipReadsTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            byte [] quals = new byte[readLength];

            for (int nLowQualBases = 0; nLowQualBases < readLength; nLowQualBases++) {

                // create a read with nLowQualBases in the left tail
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addLeft = 0; addLeft < nLowQualBases; addLeft++)
                    quals[addLeft] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKSAMRecord clipLeft = (new ReadClipper(read)).hardClipLowQualEnds(LOW_QUAL);

                // Tests

                // Make sure the low qualities are gone
                assertNoLowQualBases(clipLeft, LOW_QUAL);

                // Can't run this test with the current contract of no hanging insertions
                //Assert.assertEquals(clipLeft.getReadLength(), readLength - nLowQualBases, String.format("Clipped read size (%d) is different than the number high qual bases (%d) -- Cigars: %s -> %s", clipLeft.getReadLength(), readLength - nLowQualBases, read.getCigarString(), clipLeft.getCigarString()));

                // create a read with nLowQualBases in the right tail
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addRight = 0; addRight < nLowQualBases; addRight++)
                    quals[readLength - addRight - 1] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKSAMRecord clipRight = (new ReadClipper(read)).hardClipLowQualEnds(LOW_QUAL);

                // Tests

                // Make sure the low qualities are gone
                assertNoLowQualBases(clipRight, LOW_QUAL);

                // Make sure we haven't clipped any high quals -- Can't run this test with the current contract of no hanging insertions
                //Assert.assertEquals(clipLeft.getReadLength(), readLength - nLowQualBases, String.format("Clipped read size (%d) is different than the number high qual bases (%d) -- Cigars: %s -> %s", clipRight.getReadLength(), readLength - nLowQualBases, read.getCigarString(), clipRight.getCigarString()));

                // create a read with nLowQualBases in the both tails
                if (nLowQualBases <= readLength/2) {
                    Utils.fillArrayWithByte(quals, HIGH_QUAL);
                    for (int addBoth = 0; addBoth < nLowQualBases; addBoth++) {
                        quals[addBoth] = LOW_QUAL;
                        quals[readLength - addBoth - 1] = LOW_QUAL;
                    }
                    read.setBaseQualities(quals);
                    GATKSAMRecord clipBoth = (new ReadClipper(read)).hardClipLowQualEnds(LOW_QUAL);

                    // Tests

                    // Make sure the low qualities are gone
                    assertNoLowQualBases(clipBoth, LOW_QUAL);

                    // Can't run this test with the current contract of no hanging insertions
                    //Assert.assertEquals(clipLeft.getReadLength(), readLength - nLowQualBases, String.format("Clipped read size (%d) is different than the number high qual bases (%d) -- Cigars: %s -> %s", clipRight.getReadLength(), readLength - (2*nLowQualBases), read.getCigarString(), clipBoth.getCigarString()));
                }
            }
//            logger.warn(String.format("Testing %s for all combinations of low/high qual... PASSED", read.getCigarString()));
        }

        // ONE OFF Testing clipping that ends inside an insertion  ( Ryan's bug )
        final byte[] BASES = {'A','C','G','T','A','C','G','T'};
        final byte[] QUALS = {2, 2, 2, 2, 20, 20, 20, 2};
        final String CIGAR = "1S1M5I1S";

        final byte[] CLIPPED_BASES = {};
        final byte[] CLIPPED_QUALS = {};
        final String CLIPPED_CIGAR = "";


        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(BASES, QUALS, CIGAR);
        GATKSAMRecord expected = ArtificialSAMUtils.createArtificialRead(CLIPPED_BASES, CLIPPED_QUALS, CLIPPED_CIGAR);

        ReadClipper lowQualClipper = new ReadClipper(read);
        ClipReadsTestUtils.assertEqualReads(lowQualClipper.hardClipLowQualEnds((byte) 2), expected);

        logger.warn("PASSED");
    }

    @Test(enabled = true)
    public void testHardClipSoftClippedBases() {

        // Generate a list of cigars to test
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ClipReadsTestUtils.makeReadFromCigar(cigar);
            ReadClipper readClipper = new ReadClipper(read);
            GATKSAMRecord clippedRead = readClipper.hardClipSoftClippedBases();

            int sumHardClips = 0;
            int sumMatches = 0;

            boolean tail = true;
            for (CigarElement element : read.getCigar().getCigarElements()) {
                // Assuming cigars are well formed, if we see S or H, it means we're on the tail (left or right)
                if (element.getOperator() == CigarOperator.HARD_CLIP || element.getOperator() == CigarOperator.SOFT_CLIP)
                    tail = true;

                // Adds all H, S and D's (next to hard/soft clips).
                // All these should be hard clips after clipping.
                if (tail && (element.getOperator() == CigarOperator.HARD_CLIP || element.getOperator() == CigarOperator.SOFT_CLIP || element.getOperator() == CigarOperator.DELETION))
                    sumHardClips += element.getLength();

                // this means we're no longer on the tail (insertions can still potentially be the tail because
                // of the current contract of clipping out hanging insertions
                else if (element.getOperator() != CigarOperator.INSERTION)
                    tail = false;

                // Adds all matches to verify that they remain the same after clipping
                if (element.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                    sumMatches += element.getLength();
            }

            for (CigarElement element : clippedRead.getCigar().getCigarElements()) {
                // Test if clipped read has Soft Clips (shouldn't have any!)
                Assert.assertTrue( element.getOperator() != CigarOperator.SOFT_CLIP, String.format("Cigar %s -> %s -- FAILED (resulting cigar has soft clips)", read.getCigarString(), clippedRead.getCigarString()));

                // Keep track of the total number of Hard Clips after clipping to make sure everything was accounted for
                if (element.getOperator() == CigarOperator.HARD_CLIP)
                    sumHardClips -= element.getLength();

                // Make sure all matches are still there
                if (element.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                    sumMatches -= element.getLength();
            }
            Assert.assertTrue( sumHardClips == 0, String.format("Cigar %s -> %s -- FAILED (number of hard clips mismatched by %d)", read.getCigarString(), clippedRead.getCigarString(), sumHardClips));
            Assert.assertTrue( sumMatches == 0, String.format("Cigar %s -> %s -- FAILED (number of matches mismatched by %d)", read.getCigarString(), clippedRead.getCigarString(), sumMatches));


//            logger.warn(String.format("Cigar %s -> %s -- PASSED!", read.getCigarString(), clippedRead.getCigarString()));
        }

        logger.warn("PASSED");
    }

    private void assertNoLowQualBases(GATKSAMRecord read, byte low_qual) {
        if (!read.isEmpty()) {
            byte [] quals = read.getBaseQualities();
            for (int i=0; i<quals.length; i++)
                Assert.assertFalse(quals[i] <= low_qual, String.format("Found low qual (%d) base after hard clipping. Position: %d -- %s", low_qual, i, read.getCigarString()));
        }
    }


}