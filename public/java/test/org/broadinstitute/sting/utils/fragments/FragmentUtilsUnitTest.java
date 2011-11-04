/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.fragments;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test routines for read-backed pileup.
 */
public class FragmentUtilsUnitTest extends BaseTest {
    private static SAMFileHeader header;

    private class FragmentUtilsTest extends TestDataProvider {
        List<TestState> statesForPileup = new ArrayList<TestState>();
        List<TestState> statesForReads = new ArrayList<TestState>();

        private FragmentUtilsTest(String name, int readLen, int leftStart, int rightStart,
                                  boolean leftIsFirst, boolean leftIsNegative) {
            super(FragmentUtilsTest.class, String.format("%s-leftIsFirst:%b-leftIsNegative:%b", name, leftIsFirst, leftIsNegative));

            List<GATKSAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
            GATKSAMRecord left = pair.get(0);
            GATKSAMRecord right = pair.get(1);

            for ( int pos = leftStart; pos < rightStart + readLen; pos++) {
                boolean posCoveredByLeft = pos >= left.getAlignmentStart() && pos <= left.getAlignmentEnd();
                boolean posCoveredByRight = pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd();

                if ( posCoveredByLeft || posCoveredByRight ) {
                    List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
                    List<Integer> offsets = new ArrayList<Integer>();

                    if ( posCoveredByLeft ) {
                        reads.add(left);
                        offsets.add(pos - left.getAlignmentStart());
                    }

                    if ( posCoveredByRight ) {
                        reads.add(right);
                        offsets.add(pos - right.getAlignmentStart());
                    }

                    boolean shouldBeFragment = posCoveredByLeft && posCoveredByRight;
                    ReadBackedPileup pileup = new ReadBackedPileupImpl(null, reads, offsets);
                    TestState testState = new TestState(shouldBeFragment ? 0 : 1, shouldBeFragment ? 1 : 0, pileup, null);
                    statesForPileup.add(testState);
                }

                TestState testState = left.getAlignmentEnd() >= right.getAlignmentStart() ? new TestState(0, 1, null, pair) : new TestState(2, 0, null, pair);
                statesForReads.add(testState);
            }
        }
    }

    private class TestState {
        int expectedSingletons, expectedPairs;
        ReadBackedPileup pileup;
        List<GATKSAMRecord> rawReads;

        private TestState(final int expectedSingletons, final int expectedPairs, final ReadBackedPileup pileup, final List<GATKSAMRecord> rawReads) {
            this.expectedSingletons = expectedSingletons;
            this.expectedPairs = expectedPairs;
            this.pileup = pileup;
            this.rawReads = rawReads;
        }
    }

    @DataProvider(name = "fragmentUtilsTest")
    public Object[][] createTests() {
        for ( boolean leftIsFirst : Arrays.asList(true, false) ) {
            for ( boolean leftIsNegative : Arrays.asList(true, false) ) {
                // Overlapping pair
                // ---->        [first]
                //   <---       [second]
                new FragmentUtilsTest("overlapping-pair", 10, 1, 5, leftIsFirst, leftIsNegative);

                // Non-overlapping pair
                // ---->
                //          <----
                new FragmentUtilsTest("nonoverlapping-pair", 10, 1, 15, leftIsFirst, leftIsNegative);
            }
        }

        return FragmentUtilsTest.getTests(FragmentUtilsTest.class);
    }

    @Test(enabled = true, dataProvider = "fragmentUtilsTest")
    public void testAsPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            ReadBackedPileup rbp = testState.pileup;
            FragmentCollection<PileupElement> fp = FragmentUtils.create(rbp);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = true, dataProvider = "fragmentUtilsTest")
    public void testAsListOfReadsFromPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            FragmentCollection<GATKSAMRecord> fp = FragmentUtils.create(testState.pileup.getReads());
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = true, dataProvider = "fragmentUtilsTest")
    public void testAsListOfReads(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForReads ) {
            FragmentCollection<GATKSAMRecord> fp = FragmentUtils.create(testState.rawReads);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = true, expectedExceptions = IllegalArgumentException.class)
    public void testOutOfOrder() {
        final List<GATKSAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", 100, 1, 50, true, true);
        final GATKSAMRecord left = pair.get(0);
        final GATKSAMRecord right = pair.get(1);
        final List<GATKSAMRecord> reads = Arrays.asList(right, left); // OUT OF ORDER!
        final List<Integer> offsets = Arrays.asList(0, 50);
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(null, reads, offsets);
        FragmentUtils.create(pileup); // should throw exception
    }

    @BeforeTest
    public void setup() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
    }
}
