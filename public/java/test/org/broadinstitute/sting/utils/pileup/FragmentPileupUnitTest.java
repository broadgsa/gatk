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

package org.broadinstitute.sting.utils.pileup;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test routines for read-backed pileup.
 */
public class FragmentPileupUnitTest extends BaseTest {
    private static SAMFileHeader header;

    private class FragmentPileupTest extends TestDataProvider {
        List<TestState> states = new ArrayList<TestState>();

        private FragmentPileupTest(String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
            super(FragmentPileupTest.class, String.format("%s-leftIsFirst:%b-leftIsNegative:%b", name, leftIsFirst, leftIsNegative));

            for ( int pos = leftStart; pos < rightStart + readLen; pos++) {
                SAMRecord left = ArtificialSAMUtils.createArtificialRead(header, "readpair", 0, leftStart, readLen);
                SAMRecord right = ArtificialSAMUtils.createArtificialRead(header, "readpair", 0, rightStart, readLen);

                left.setProperPairFlag(true);
                right.setProperPairFlag(true);

                left.setFirstOfPairFlag(leftIsFirst);
                right.setFirstOfPairFlag(! leftIsFirst);

                left.setReadNegativeStrandFlag(leftIsNegative);
                left.setMateNegativeStrandFlag(!leftIsNegative);
                right.setReadNegativeStrandFlag(!leftIsNegative);
                right.setMateNegativeStrandFlag(leftIsNegative);

                left.setMateAlignmentStart(right.getAlignmentStart());
                right.setMateAlignmentStart(left.getAlignmentStart());

                left.setMateReferenceIndex(0);
                right.setMateReferenceIndex(0);

                int isize = rightStart + readLen - leftStart;
                left.setInferredInsertSize(isize);
                right.setInferredInsertSize(-isize);

                boolean posCoveredByLeft = pos >= left.getAlignmentStart() && pos <= left.getAlignmentEnd();
                boolean posCoveredByRight = pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd();

                if ( posCoveredByLeft || posCoveredByRight ) {
                    List<SAMRecord> reads = new ArrayList<SAMRecord>();
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
                    TestState testState = new TestState(shouldBeFragment, pileup);
                    states.add(testState);
                }
            }
        }
    }

    private class TestState {
        boolean shouldBeFragment;
        ReadBackedPileup pileup;

        private TestState(final boolean shouldBeFragment, final ReadBackedPileup pileup) {
            this.shouldBeFragment = shouldBeFragment;
            this.pileup = pileup;
        }
    }

    @DataProvider(name = "fragmentPileupTest")
    public Object[][] createTests() {
        for ( boolean leftIsFirst : Arrays.asList(true, false) ) {
            for ( boolean leftIsNegative : Arrays.asList(true, false) ) {
                // Overlapping pair
                // ---->        [first]
                //   <---       [second]
                new FragmentPileupTest("overlapping-pair", 10, 1, 5, leftIsFirst, leftIsNegative);

                // Non-overlapping pair
                // ---->
                //          <----
                new FragmentPileupTest("nonoverlapping-pair", 10, 1, 15, leftIsFirst, leftIsNegative);
            }
        }

        return FragmentPileupTest.getTests(FragmentPileupTest.class);
    }

    @Test(enabled = true, dataProvider = "fragmentPileupTest")
    public void testMe(FragmentPileupTest test) {
        for ( TestState testState : test.states ) {
            ReadBackedPileup rbp = testState.pileup;
            FragmentPileup fp = new FragmentPileup(rbp);
            Assert.assertEquals(fp.getTwoReadPileup().size(), testState.shouldBeFragment ? 1 : 0);
            Assert.assertEquals(fp.getOneReadPileup().size(), testState.shouldBeFragment ? 0 : 1);
        }
    }


    @BeforeTest
    public void setup() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
    }
}
