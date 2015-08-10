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

package org.broadinstitute.gatk.engine.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.downsampling.FractionalDownsampler;
import org.broadinstitute.gatk.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class FractionalDownsamplerUnitTest extends BaseTest {

    private static class FractionalDownsamplerTest extends TestDataProvider {
        double fraction;
        int totalReads;
        int expectedMinNumReadsAfterDownsampling;
        int expectedMaxNumReadsAfterDownsampling;
        int expectedMinDiscardedItems;
        int expectedMaxDiscardedItems;

        private static final double EXPECTED_ACCURACY = 0.05; // should be accurate to within +/- this percent

        public FractionalDownsamplerTest( double fraction, int totalReads ) {
            super(FractionalDownsamplerTest.class);

            this.fraction = fraction;
            this.totalReads = totalReads;

            calculateExpectations();

            setName(String.format("%s: fraction=%.2f totalReads=%d expectedMinNumReadsAfterDownsampling=%d expectedMaxNumReadsAfterDownsampling=%d",
                    getClass().getSimpleName(), fraction, totalReads, expectedMinNumReadsAfterDownsampling, expectedMaxNumReadsAfterDownsampling));
        }

        private void calculateExpectations() {
            // Require an exact match in the 0% and 100% cases
            if ( fraction == 0.0 ) {
                expectedMinNumReadsAfterDownsampling = expectedMaxNumReadsAfterDownsampling = 0;
                expectedMinDiscardedItems = expectedMaxDiscardedItems = totalReads;
            }
            else if ( fraction == 1.0 ) {
                expectedMinNumReadsAfterDownsampling = expectedMaxNumReadsAfterDownsampling = totalReads;
                expectedMinDiscardedItems = expectedMaxDiscardedItems = 0;
            }
            else {
                expectedMinNumReadsAfterDownsampling = Math.max((int)((fraction - EXPECTED_ACCURACY) * totalReads), 0);
                expectedMaxNumReadsAfterDownsampling = Math.min((int) ((fraction + EXPECTED_ACCURACY) * totalReads), totalReads);
                expectedMinDiscardedItems = totalReads - expectedMaxNumReadsAfterDownsampling;
                expectedMaxDiscardedItems = totalReads - expectedMinNumReadsAfterDownsampling;
            }
        }

        public Collection<SAMRecord> createReads() {
            Collection<SAMRecord> reads = new ArrayList<SAMRecord>(totalReads);

            SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
            reads.addAll(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(totalReads, header, "foo", 0, 1, 100));

            return reads;
        }
    }

    @DataProvider(name = "FractionalDownsamplerTestDataProvider")
    public Object[][] createFractionalDownsamplerTestData() {
        for ( double fraction : Arrays.asList(0.0, 0.25, 0.5, 0.75, 1.0) ) {
            for ( int totalReads : Arrays.asList(0, 1000, 10000) ) {
                new FractionalDownsamplerTest(fraction, totalReads);
            }
        }

        return FractionalDownsamplerTest.getTests(FractionalDownsamplerTest.class);
    }

    @Test(dataProvider = "FractionalDownsamplerTestDataProvider")
    public void runFractionalDownsamplerTest( FractionalDownsamplerTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();

        ReadsDownsampler<SAMRecord> downsampler = new FractionalDownsampler<SAMRecord>(test.fraction);

        downsampler.submit(test.createReads());

        if ( test.totalReads > 0 ) {
            if ( test.fraction > FractionalDownsamplerTest.EXPECTED_ACCURACY ) {
                Assert.assertTrue(downsampler.hasFinalizedItems());
                Assert.assertTrue(downsampler.peekFinalized() != null);
            }
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        downsampler.signalEndOfInput();

        if ( test.totalReads > 0 ) {
            if ( test.fraction > FractionalDownsamplerTest.EXPECTED_ACCURACY ) {
                Assert.assertTrue(downsampler.hasFinalizedItems());
                Assert.assertTrue(downsampler.peekFinalized() != null);
            }
            Assert.assertFalse(downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekPending() == null);
        }
        else {
            Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
            Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);
        }

        List<SAMRecord> downsampledReads = downsampler.consumeFinalizedItems();
        Assert.assertFalse(downsampler.hasFinalizedItems() || downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.peekFinalized() == null && downsampler.peekPending() == null);

        Assert.assertTrue(downsampledReads.size() >= test.expectedMinNumReadsAfterDownsampling &&
                          downsampledReads.size() <= test.expectedMaxNumReadsAfterDownsampling);

        Assert.assertTrue(downsampler.getNumberOfDiscardedItems() >= test.expectedMinDiscardedItems &&
                          downsampler.getNumberOfDiscardedItems() <= test.expectedMaxDiscardedItems);

        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), test.totalReads - downsampledReads.size());

        downsampler.resetStats();
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }
}
