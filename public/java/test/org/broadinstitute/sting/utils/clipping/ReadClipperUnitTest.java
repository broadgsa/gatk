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

package org.broadinstitute.sting.utils.clipping;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.List;

/**
 * User: roger
 * Date: 9/28/11
 */
public class ReadClipperUnitTest extends BaseTest {

    List<Cigar> cigarList;
    int maximumCigarSize = 6;                                                                                           // 6 is the minimum necessary number to try all combinations of cigar types with guarantee of clipping an element with length = 2

    @BeforeClass
    public void init() {
        cigarList = ReadClipperTestUtils.generateCigarList(maximumCigarSize);
    }

    @Test(enabled = true)
    public void testHardClipBothEndsByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                GATKSAMRecord clippedRead = ReadClipper.hardClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getAlignmentStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
                Assert.assertTrue(clippedRead.getAlignmentEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
                assertUnclippedLimits(read, clippedRead);
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipByReadCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            for (int i = 0; i < readLength; i++) {
                GATKSAMRecord clipLeft = ReadClipper.hardClipByReadCoordinates(read, 0, i);
                Assert.assertTrue(clipLeft.getReadLength() <= readLength - i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipLeft.getCigarString()));
                assertUnclippedLimits(read, clipLeft);

                GATKSAMRecord clipRight = ReadClipper.hardClipByReadCoordinates(read, i, readLength - 1);
                Assert.assertTrue(clipRight.getReadLength() <= i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipRight.getCigarString()));
                assertUnclippedLimits(read, clipRight);
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int start = read.getSoftStart();
            int stop = read.getSoftEnd();

            for (int i = start; i <= stop; i++) {
                GATKSAMRecord clipLeft = (new ReadClipper(read)).hardClipByReferenceCoordinates(-1, i);
                if (!clipLeft.isEmpty()) {
                    Assert.assertTrue(clipLeft.getAlignmentStart() >= Math.min(read.getAlignmentEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getAlignmentStart(), i + 1, read.getCigarString(), clipLeft.getCigarString()));
                    assertUnclippedLimits(read, clipLeft);
                }

                GATKSAMRecord clipRight = (new ReadClipper(read)).hardClipByReferenceCoordinates(i, -1);
                if (!clipRight.isEmpty() && clipRight.getAlignmentStart() <= clipRight.getAlignmentEnd()) {             // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                    Assert.assertTrue(clipRight.getAlignmentEnd() <= Math.max(read.getAlignmentStart(), i - 1), String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getAlignmentEnd(), i - 1, read.getCigarString(), clipRight.getCigarString()));
                    assertUnclippedLimits(read, clipRight);
                }
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesLeftTail() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            if (read.getSoftStart() == alnStart) {                                                                      // we can't test left clipping if the read has hanging soft clips on the left side
                for (int i = alnStart; i <= alnEnd; i++) {
                    GATKSAMRecord clipLeft = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, i);

                    if (!clipLeft.isEmpty()) {
//                        System.out.println(String.format("Left Tail [%d]: %s (%d,%d,%d : %d,%d,%d) -> %s (%d,%d,%d : %d,%d,%d)", i, cigar.toString(), read.getUnclippedStart(), read.getSoftStart(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getSoftEnd(), read.getUnclippedEnd(), clipLeft.getCigarString(), clipLeft.getUnclippedStart(), clipLeft.getSoftStart(), clipLeft.getAlignmentStart(), clipLeft.getAlignmentEnd(), clipLeft.getSoftEnd(), clipLeft.getUnclippedEnd()));
                        Assert.assertTrue(clipLeft.getAlignmentStart() >= i + 1, String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getAlignmentStart(), i + 1, read.getCigarString(), clipLeft.getCigarString()));
                        assertUnclippedLimits(read, clipLeft);
                    }
                }
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesRightTail() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            if (read.getSoftEnd() == alnEnd) {                                                                          // we can't test right clipping if the read has hanging soft clips on the right side
                for (int i = alnStart; i <= alnEnd; i++) {
                    GATKSAMRecord clipRight = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, i);
                    if (!clipRight.isEmpty() && clipRight.getAlignmentStart() <= clipRight.getAlignmentEnd()) {  // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                        Assert.assertTrue(clipRight.getAlignmentEnd() <= i - 1, String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getAlignmentEnd(), i - 1, read.getCigarString(), clipRight.getCigarString()));
                        assertUnclippedLimits(read, clipRight);
                    }
                }
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipLowQualEnds() {
        final byte LOW_QUAL = 2;
        final byte HIGH_QUAL = 30;

        /** create a read for every cigar permutation */
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            byte[] quals = new byte[readLength];

            for (int nLowQualBases = 0; nLowQualBases < readLength; nLowQualBases++) {

                /**  create a read with nLowQualBases in the left tail */
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addLeft = 0; addLeft < nLowQualBases; addLeft++)
                    quals[addLeft] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKSAMRecord clipLeft = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipLeft, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases in the right tail */
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addRight = 0; addRight < nLowQualBases; addRight++)
                    quals[readLength - addRight - 1] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKSAMRecord clipRight = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipRight, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases on both tails */
                if (nLowQualBases <= readLength / 2) {
                    Utils.fillArrayWithByte(quals, HIGH_QUAL);
                    for (int addBoth = 0; addBoth < nLowQualBases; addBoth++) {
                        quals[addBoth] = LOW_QUAL;
                        quals[readLength - addBoth - 1] = LOW_QUAL;
                    }
                    read.setBaseQualities(quals);
                    GATKSAMRecord clipBoth = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                    checkClippedReadsForLowQualEnds(read, clipBoth, LOW_QUAL, 2*nLowQualBases);
                }
            }
        }
    }

    @Test(enabled = true)
    public void testHardClipSoftClippedBases() {
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            GATKSAMRecord clippedRead = ReadClipper.hardClipSoftClippedBases(read);
            CigarCounter original = new CigarCounter(read);
            CigarCounter clipped = new CigarCounter(clippedRead);

            assertUnclippedLimits(read, clippedRead);                                                                   // Make sure limits haven't changed
            original.assertHardClippingSoftClips(clipped);                                                              // Make sure we have only clipped SOFT_CLIPS
        }
    }

    @Test(enabled = false)
    public void testHardClipLeadingInsertions() {
        for (Cigar cigar : cigarList) {
            if (startsWithInsertion(cigar)) {
                GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
                GATKSAMRecord clippedRead = ReadClipper.hardClipLeadingInsertions(read);

                assertUnclippedLimits(read, clippedRead);        // Make sure limits haven't changed

                int expectedLength = read.getReadLength() - leadingCigarElementLength(read.getCigar(), CigarOperator.INSERTION);
                if (cigarHasElementsDifferentThanInsertionsAndHardClips(read.getCigar()))
                    expectedLength -= leadingCigarElementLength(ReadClipperTestUtils.invertCigar(read.getCigar()), CigarOperator.INSERTION);

                if (!clippedRead.isEmpty()) {
                    Assert.assertEquals(expectedLength, clippedRead.getReadLength(), String.format("%s -> %s", read.getCigarString(), clippedRead.getCigarString()));  // check that everything else is still there
                    Assert.assertFalse(startsWithInsertion(clippedRead.getCigar()));                                                                                   // check that the insertions are gone
                } else
                    Assert.assertTrue(expectedLength == 0, String.format("expected length: %d", expectedLength));                                                      // check that the read was expected to be fully clipped
            }
        }
    }

    @Test(enabled = true)
    public void testRevertSoftClippedBases() {
        for (Cigar cigar : cigarList) {
            final int leadingSoftClips = leadingCigarElementLength(cigar, CigarOperator.SOFT_CLIP);
            final int tailSoftClips = leadingCigarElementLength(ReadClipperTestUtils.invertCigar(cigar), CigarOperator.SOFT_CLIP);

            final GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final GATKSAMRecord unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed

            if (leadingSoftClips > 0 || tailSoftClips > 0) {
                final int expectedStart = read.getAlignmentStart() - leadingSoftClips;
                final int expectedEnd = read.getAlignmentEnd() + tailSoftClips;

                Assert.assertEquals(unclipped.getAlignmentStart(), expectedStart);
                Assert.assertEquals(unclipped.getAlignmentEnd(), expectedEnd);
            } else
                Assert.assertEquals(read.getCigarString(), unclipped.getCigarString());
        }
    }

    @Test(enabled = true)
    public void testRevertSoftClippedBasesWithThreshold() {
        for (Cigar cigar : cigarList) {
            final int leadingSoftClips = leadingCigarElementLength(cigar, CigarOperator.SOFT_CLIP);
            final int tailSoftClips = leadingCigarElementLength(ReadClipperTestUtils.invertCigar(cigar), CigarOperator.SOFT_CLIP);

            final GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final GATKSAMRecord unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed
            Assert.assertNull(read.getCigar().isValid(null, -1));
            Assert.assertNull(unclipped.getCigar().isValid(null, -1));

            if (!(leadingSoftClips > 0 || tailSoftClips > 0))
                Assert.assertEquals(read.getCigarString(), unclipped.getCigarString());

        }
    }


    private void assertNoLowQualBases(GATKSAMRecord read, byte low_qual) {
        if (!read.isEmpty()) {
            byte[] quals = read.getBaseQualities();
            for (int i = 0; i < quals.length; i++)
                Assert.assertFalse(quals[i] <= low_qual, String.format("Found low qual (%d) base after hard clipping. Position: %d -- %s", low_qual, i, read.getCigarString()));
        }
    }

    private void checkClippedReadsForLowQualEnds(GATKSAMRecord read, GATKSAMRecord clippedRead, byte lowQual, int nLowQualBases) {
        assertUnclippedLimits(read, clippedRead);                                                                       // Make sure limits haven't changed
        assertNoLowQualBases(clippedRead, lowQual);                                                                     // Make sure the low qualities are gone
        assertNumberOfBases(read, clippedRead, nLowQualBases);                                                          // Make sure only low quality bases were clipped
    }

    /**
     * Asserts that clipping doesn't change the getUnclippedStart / getUnclippedEnd
     *
     * @param original original read
     * @param clipped clipped read
     */
    private void assertUnclippedLimits(GATKSAMRecord original, GATKSAMRecord clipped) {
        if (ReadClipperTestUtils.readHasNonClippedBases(clipped)) {
            Assert.assertEquals(original.getUnclippedStart(), clipped.getUnclippedStart());
            Assert.assertEquals(original.getUnclippedEnd(), clipped.getUnclippedEnd());
        }
    }

    private void assertNumberOfBases(GATKSAMRecord read, GATKSAMRecord clipLeft, int nLowQualBases) {
        if (read.getCigarString().contains("M"))
            Assert.assertEquals(clipLeft.getReadLength(), read.getReadLength() - nLowQualBases, String.format("Clipped read size (%d) is different than the number high qual bases (%d) -- Cigars: %s -> %s", clipLeft.getReadLength(), read.getReadLength() - nLowQualBases, read.getCigarString(), clipLeft.getCigarString()));
    }


    private boolean startsWithInsertion(Cigar cigar) {
        return leadingCigarElementLength(cigar, CigarOperator.INSERTION) > 0;
    }

    private int leadingCigarElementLength(Cigar cigar, CigarOperator operator) {
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator() == operator)
                return cigarElement.getLength();
            if (cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                break;
        }
        return 0;
    }

    private boolean cigarHasElementsDifferentThanInsertionsAndHardClips(Cigar cigar) {
        for (CigarElement cigarElement : cigar.getCigarElements())
            if (cigarElement.getOperator() != CigarOperator.INSERTION && cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                return true;
        return false;
    }

    private class CigarCounter {
        private HashMap<CigarOperator, Integer> counter;

        public Integer getCounterForOp(CigarOperator operator) {
            return counter.get(operator);
        }

        public CigarCounter(GATKSAMRecord read) {
            CigarOperator[] operators = CigarOperator.values();
            counter = new HashMap<CigarOperator, Integer>(operators.length);

            for (CigarOperator op : operators)
                counter.put(op, 0);

            for (CigarElement cigarElement : read.getCigar().getCigarElements())
                counter.put(cigarElement.getOperator(), counter.get(cigarElement.getOperator()) + cigarElement.getLength());
        }

        public boolean assertHardClippingSoftClips(CigarCounter clipped) {
            for (CigarOperator op : counter.keySet()) {
                if (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP) {
                    int counterTotal = counter.get(CigarOperator.HARD_CLIP) + counter.get(CigarOperator.SOFT_CLIP);
                    int clippedHard = clipped.getCounterForOp(CigarOperator.HARD_CLIP);
                    int clippedSoft = clipped.getCounterForOp(CigarOperator.SOFT_CLIP);

                    Assert.assertEquals(counterTotal, clippedHard);
                    Assert.assertTrue(clippedSoft == 0);
                } else
                    Assert.assertEquals(counter.get(op), clipped.getCounterForOp(op));
            }
            return true;
        }

    }

}