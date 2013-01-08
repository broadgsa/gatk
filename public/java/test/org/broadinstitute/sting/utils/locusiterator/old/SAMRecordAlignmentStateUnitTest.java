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

package org.broadinstitute.sting.utils.locusiterator.old;

import org.broadinstitute.sting.utils.locusiterator.LIBS_position;
import org.broadinstitute.sting.utils.locusiterator.LocusIteratorByStateBaseTest;
import org.broadinstitute.sting.utils.locusiterator.old.SAMRecordAlignmentState;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class SAMRecordAlignmentStateUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "AlignmentStateTest")
    public Object[][] makeAlignmentStateTest() {
//        return new Object[][]{{new LIBSTest("1I", 1)}};
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "AlignmentStateTest")
    public void testAlignmentStateTest(LIBSTest params) {
        final GATKSAMRecord read = params.makeRead();
        final SAMRecordAlignmentState state = new SAMRecordAlignmentState(read);
        final LIBS_position tester = new LIBS_position(read);

        Assert.assertSame(state.getRead(), read);
        Assert.assertNotNull(state.toString());

        int bpVisited = 0;
        int lastOffset = -1;
        while ( state.stepForwardOnGenome() != null ) {
            bpVisited++;
            tester.stepForwardOnGenome();
            Assert.assertTrue(state.getReadOffset() >= lastOffset, "Somehow read offsets are decreasing: lastOffset " + lastOffset + " current " + state.getReadOffset());
            Assert.assertEquals(state.getReadOffset(), tester.getCurrentReadOffset(), "Read offsets are wrong at " + bpVisited);

            // TODO -- state.peekBackwardOnGenome();
            // TODO -- state.peekForwardOnGenome();
            // TODO -- state.getCurrentCigarOperator()
            // TODO -- state.getGenomeOffset();
            // TODO -- state.getGenomePosition();
            // TODO -- Assert.assertEquals(state.getLocation(genomeLocParser), EXPECTATION);

            lastOffset = state.getReadOffset();
        }

        // min is one because always visit something, even for 10I reads
        final int expectedBpToVisit = read.getAlignmentEnd() - read.getAlignmentStart() + 1;
        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
    }
}
