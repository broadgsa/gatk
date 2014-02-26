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

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.util.*;

public class LevelingDownsamplerUnitTest extends BaseTest {

    private static class LevelingDownsamplerUniformStacksTest extends TestDataProvider {
        public enum DataStructure { LINKED_LIST, ARRAY_LIST }

        int targetSize;
        int numStacks;
        int stackSize;
        DataStructure dataStructure;
        int expectedSize;

        public LevelingDownsamplerUniformStacksTest( int targetSize, int numStacks, int stackSize, DataStructure dataStructure ) {
            super(LevelingDownsamplerUniformStacksTest.class);

            this.targetSize = targetSize;
            this.numStacks = numStacks;
            this.stackSize = stackSize;
            this.dataStructure = dataStructure;
            expectedSize = calculateExpectedDownsampledStackSize();

            setName(String.format("%s: targetSize=%d numStacks=%d stackSize=%d dataStructure=%s expectedSize=%d",
                    getClass().getSimpleName(), targetSize, numStacks, stackSize, dataStructure, expectedSize));
        }

        public Collection<List<Object>> createStacks() {
            Collection<List<Object>> stacks = new ArrayList<List<Object>>();

            for ( int i = 1; i <= numStacks; i++ ) {
                List<Object> stack = dataStructure == DataStructure.LINKED_LIST ? new LinkedList<Object>() : new ArrayList<Object>();

                for ( int j = 1; j <= stackSize; j++ ) {
                    stack.add(new Object());
                }

                stacks.add(stack);
            }

            return stacks;
        }

        private int calculateExpectedDownsampledStackSize() {
            int numItemsToRemove = numStacks * stackSize - targetSize;

            if ( numStacks == 0 ) {
                return 0;
            }
            else if ( numItemsToRemove <= 0 ) {
                return stackSize;
            }

            return Math.max(1, stackSize - (numItemsToRemove / numStacks));
        }
    }

    @DataProvider(name = "UniformStacksDataProvider")
    public Object[][] createUniformStacksTestData() {
        for ( int targetSize = 1; targetSize <= 10000; targetSize *= 10 ) {
            for ( int numStacks = 0; numStacks <= 10; numStacks++ ) {
                for ( int stackSize = 1; stackSize <= 1000; stackSize *= 10 ) {
                    for ( LevelingDownsamplerUniformStacksTest.DataStructure dataStructure : LevelingDownsamplerUniformStacksTest.DataStructure.values() ) {
                        new LevelingDownsamplerUniformStacksTest(targetSize, numStacks, stackSize, dataStructure);
                    }
                }
            }
        }

        return LevelingDownsamplerUniformStacksTest.getTests(LevelingDownsamplerUniformStacksTest.class);
    }

    @Test( dataProvider = "UniformStacksDataProvider" )
    public void testLevelingDownsamplerWithUniformStacks( LevelingDownsamplerUniformStacksTest test ) {
        logger.warn("Running test: " + test);

        GenomeAnalysisEngine.resetRandomGenerator();

        Downsampler<List<Object>> downsampler = new LevelingDownsampler<List<Object>, Object>(test.targetSize);

        downsampler.submit(test.createStacks());

        if ( test.numStacks > 0 ) {
            Assert.assertFalse(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() == null);
            Assert.assertTrue(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() != null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.numStacks > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        final int sizeFromDownsampler = downsampler.size();
        List<List<Object>> downsampledStacks = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertEquals(downsampledStacks.size(), test.numStacks);

        int totalRemainingItems = 0;
        for ( List<Object> stack : downsampledStacks ) {
            Assert.assertTrue(Math.abs(stack.size() - test.expectedSize) <= 1);
            totalRemainingItems += stack.size();
        }

        Assert.assertEquals(sizeFromDownsampler, totalRemainingItems);
        int numItemsReportedDiscarded = downsampler.getNumberOfDiscardedItems();
        int numItemsActuallyDiscarded = test.numStacks * test.stackSize - totalRemainingItems;

        Assert.assertEquals(numItemsReportedDiscarded, numItemsActuallyDiscarded);

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);

        Assert.assertTrue(totalRemainingItems <= Math.max(test.targetSize, test.numStacks));
    }

    @Test
    public void testDoNotDiscardReducedReads() {
        GenomeAnalysisEngine.resetRandomGenerator();
        final Downsampler<LinkedList<AlignmentStateMachine>> downsampler = new LevelingDownsampler<LinkedList<AlignmentStateMachine>, AlignmentStateMachine>(1);

        final Collection<LinkedList<AlignmentStateMachine>> groups = new LinkedList<LinkedList<AlignmentStateMachine>>();
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
        final int[] baseCounts = { 10, 10, 10, 10, 10 };

        for ( int alignmentStart : Arrays.asList(1, 2, 3) ) {
            final LinkedList<AlignmentStateMachine> group = new LinkedList<AlignmentStateMachine>();
            for ( int i = 1; i <= 10; i++ ) {
                group.add(new AlignmentStateMachine(ArtificialSAMUtils.createArtificialReducedRead(header, "foo", 0, alignmentStart, 5, baseCounts)));
            }
            groups.add(group);
        }

        downsampler.submit(groups);
        downsampler.signalEndOfInput();

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0, "wrong number of items discarded by the downsampler");
        Assert.assertTrue(downsampler.hasFinalizedItems(), "downsampler should have finalized items but doesn't");
        Assert.assertEquals(downsampler.size(), 30, "downsampler size() reports wrong number of items");

        final Collection<LinkedList<AlignmentStateMachine>> groupsReturned = downsampler.consumeFinalizedItems();

        Assert.assertEquals(groupsReturned.size(), 3, "wrong number of groups returned by the downsampler");

        for ( LinkedList<AlignmentStateMachine> group : groupsReturned ) {
            Assert.assertEquals(group.size(), 10, "group has wrong size after downsampling");

            for ( AlignmentStateMachine state : group ) {
                Assert.assertTrue(state.isReducedRead());
            }
        }
    }
}
