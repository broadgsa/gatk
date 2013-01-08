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

package org.broadinstitute.sting.utils.locusiterator;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class AlignmentStateMachineUnitTest extends LocusIteratorByStateBaseTest {
    @DataProvider(name = "AlignmentStateMachineTest")
    public Object[][] makeAlignmentStateMachineTest() {
//        return new Object[][]{{new LIBSTest("2X2D2P2X", 1)}};
//        return createLIBSTests(
//                Arrays.asList(1, 2),
//                Arrays.asList(5));
        return createLIBSTests(
                Arrays.asList(1, 2),
                Arrays.asList(1, 2, 3, 4));
    }

    @Test(dataProvider = "AlignmentStateMachineTest")
    public void testAlignmentStateMachineTest(LIBSTest params) {
        final GATKSAMRecord read = params.makeRead();
        final AlignmentStateMachine stateMachine = new AlignmentStateMachine(read);
        final LIBS_position tester = new LIBS_position(read);

        // min is one because always visit something, even for 10I reads
        final int expectedBpToVisit = read.getAlignmentEnd() - read.getAlignmentStart() + 1;

        Assert.assertSame(stateMachine.getRead(), read);
        Assert.assertNotNull(stateMachine.toString());

        int bpVisited = 0;
        int lastOffset = -1;

        // TODO -- test state machine state before first step?

        while ( stateMachine.stepForwardOnGenome() != null ) {
            tester.stepForwardOnGenome();
            final AlignmentState state = stateMachine.getCurrent();

            Assert.assertTrue(state.getReadOffset() >= lastOffset, "Somehow read offsets are decreasing: lastOffset " + lastOffset + " current " + state.getReadOffset());
            Assert.assertEquals(state.getReadOffset(), tester.getCurrentReadOffset(), "Read offsets are wrong at " + bpVisited);

            if ( bpVisited == 0 ) {
                Assert.assertTrue(state.getPrev().isEdge());
                Assert.assertTrue(state.prevIsEdge());
            }

            if ( bpVisited == expectedBpToVisit ) {
                Assert.assertTrue(state.hasNext());
                Assert.assertTrue(state.nextIsEdge());
            }

            if ( ! state.nextIsEdge() )
                Assert.assertSame(state.getNext().getPrev(), state);

            testSequencialStatesAreConsistent(state.getPrev(), state);
            testSequencialStatesAreConsistent(state, state.getNext());

            if ( ! workAroundOpsBetweenDeletion(state.getBetweenPrevPosition()))
                Assert.assertEquals(state.isAfterDeletion(), tester.isAfterDeletedBase, "fails after deletion");
            if ( ! workAroundOpsBetweenDeletion(state.getBetweenNextPosition()))
                Assert.assertEquals(state.isBeforeDeletion(), tester.isBeforeDeletedBase, "fails before deletion");
            Assert.assertEquals(state.isAfterInsertion(), tester.isAfterInsertion, "fails after insertion");
            Assert.assertEquals(state.isBeforeInsertion(), tester.isBeforeInsertion, "Fails before insertion");
            Assert.assertEquals(state.isNextToSoftClip(), tester.isNextToSoftClip, "Fails soft clip test");

            // TODO -- fixme
            //Assert.assertEquals(state.getCigarElementCounter(), tester.currentOperatorIndex, "CigarElement indice failure");

            // TODO -- state.getGenomeOffset();
            // TODO -- state.getGenomePosition();
            // TODO -- Assert.assertEquals(state.getLocation(genomeLocParser), EXPECTATION);

            lastOffset = state.getReadOffset();
            bpVisited++;
        }

        Assert.assertTrue(stateMachine.getCurrent().isEdge());
        Assert.assertFalse(stateMachine.getCurrent().hasNext());
        Assert.assertEquals(stateMachine.getCurrent().getNext(), null);

        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
    }

    /**
     * Work around inadequate tests that aren't worth fixing.
     *
     * Look at the CIGAR 2M2P2D2P2M.  Both M states border a deletion, separated by P (padding elements).  So
     * the right answer for deletions here is true for isBeforeDeletion() and isAfterDeletion() for the first
     * and second M.  But the LIBS_position doesn't say so.
     *
     * @param elements
     * @return
     */
    private boolean workAroundOpsBetweenDeletion(final List<CigarElement> elements) {
        for ( final CigarElement elt : elements )
            if ( elt.getOperator() == CigarOperator.P || elt.getOperator() == CigarOperator.H || elt.getOperator() == CigarOperator.S )
                return true;
        return false;
    }

    private void testSequencialStatesAreConsistent(final AlignmentState left, final AlignmentState right) {
        Assert.assertSame(left.getNext(), right);
        Assert.assertSame(right.getPrev(), left);
        Assert.assertSame(left.getBetweenNextPosition(), right.getBetweenPrevPosition());
    }
}
