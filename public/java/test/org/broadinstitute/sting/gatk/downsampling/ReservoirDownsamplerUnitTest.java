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
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ReservoirDownsamplerUnitTest extends BaseTest {

    private static class ReservoirDownsamplerTest extends TestDataProvider {
        int reservoirSize;
        int totalReads;
        int expectedNumReadsAfterDownsampling;
        int expectedNumDiscardedItems;

        public ReservoirDownsamplerTest( int reservoirSize, int totalReads ) {
            super(ReservoirDownsamplerTest.class);

            this.reservoirSize = reservoirSize;
            this.totalReads = totalReads;

            expectedNumReadsAfterDownsampling = Math.min(reservoirSize, totalReads);
            expectedNumDiscardedItems = totalReads <= reservoirSize ? 0 : totalReads - reservoirSize;

            setName(String.format("%s: reservoirSize=%d totalReads=%d expectedNumReadsAfterDownsampling=%d expectedNumDiscardedItems=%d",
                    getClass().getSimpleName(), reservoirSize, totalReads, expectedNumReadsAfterDownsampling, expectedNumDiscardedItems));
        }

        public Collection<SAMRecord> createReads() {
            Collection<SAMRecord> reads = new ArrayList<SAMRecord>(totalReads);

            SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
            reads.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(totalReads, header, "foo", 0, 1, 100));

            return reads;
        }
    }

    @DataProvider(name = "ReservoirDownsamplerTestDataProvider")
    public Object[][] createReservoirDownsamplerTestData() {
        for ( int reservoirSize = 1; reservoirSize <= 10000; reservoirSize *= 10 ) {
            new ReservoirDownsamplerTest(reservoirSize, 0);
            for ( int totalReads = 1; totalReads <= 10000; totalReads *= 10 ) {
                new ReservoirDownsamplerTest(reservoirSize, totalReads);
            }
        }

        return ReservoirDownsamplerTest.getTests(ReservoirDownsamplerTest.class);
    }

    @Test(dataProvider = "ReservoirDownsamplerTestDataProvider")
    public void testReservoirDownsampler( ReservoirDownsamplerTest test ) {
        logger.warn("Running test: " + test);

        GenomeAnalysisEngine.resetRandomGenerator();

        ReadsDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(test.reservoirSize);

        downsampler.submit(test.createReads());

        if ( test.totalReads > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.totalReads > 0 ) {
            Assert.assertTrue(downsampler.hasFinalizedItems());
            Assert.assertTrue(downsampler.peekFinalized() != null);
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        Assert.assertEquals(downsampler.size(), test.expectedNumReadsAfterDownsampling);
        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertEquals(downsampledReads.size(), test.expectedNumReadsAfterDownsampling);

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), test.expectedNumDiscardedItems);
        Assert.assertEquals(test.totalReads - downsampledReads.size(), test.expectedNumDiscardedItems);

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }

    @Test
    public void testDoNotDiscardReducedReads() {
        GenomeAnalysisEngine.resetRandomGenerator();
        final ReadsDownsampler<GATKSAMRecord> downsampler = new ReservoirDownsampler<GATKSAMRecord>(1);

        final Collection<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
        final int[] baseCounts = { 10, 10, 10, 10, 10 };

        for ( int i = 1; i <= 10; i++ ) {
            reads.add(ArtificialSAMUtils.createArtificialReducedRead(header, "foo", 0, 1, 5, baseCounts));
        }
        for ( int i = 1; i <= 5; i++ ) {
            reads.add(ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 1, 5));
        }

        downsampler.submit(reads);
        downsampler.signalEndOfInput();

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 4, "wrong number of items discarded by the downsampler");
        Assert.assertTrue(downsampler.hasFinalizedItems(), "downsampler should have finalized items but doesn't");
        Assert.assertEquals(downsampler.size(), 11, "downsampler size() reports wrong number of items");

        final Collection<GATKSAMRecord> readsReturned = downsampler.consumeFinalizedItems();

        Assert.assertEquals(readsReturned.size(), 11, "wrong number of items returned by the downsampler");

        int numReducedReadsReturned = 0;
        int numNormalReadsReturned = 0;
        for ( GATKSAMRecord readReturned : readsReturned ) {
            if ( readReturned.isReducedRead() ) {
                numReducedReadsReturned++;
            }
            else {
                numNormalReadsReturned++;
            }
        }

        Assert.assertEquals(numReducedReadsReturned, 10, "wrong number of reduced reads returned by the downsampler");
        Assert.assertEquals(numNormalReadsReturned, 1, "wrong number of non-reduced reads returned by the downsampler");
    }
}
