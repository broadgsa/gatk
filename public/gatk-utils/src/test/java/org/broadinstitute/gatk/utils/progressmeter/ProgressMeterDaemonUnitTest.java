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

package org.broadinstitute.gatk.utils.progressmeter;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * UnitTests for the ProgressMeterDaemon
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ProgressMeterDaemonUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        genomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(new File(b37KGReference)));
    }

    // capture and count calls to progress
    private class TestingProgressMeter extends ProgressMeter {
        final List<Long> progressCalls = new LinkedList<Long>();

        private TestingProgressMeter(final long poll) {
            super(null, "test", new GenomeLocSortedSet(genomeLocParser), poll);
            super.start();
        }

        @Override
        protected synchronized void printProgress(boolean mustPrint) {
            progressCalls.add(System.currentTimeMillis());
        }
    }

    @DataProvider(name = "PollingData")
    public Object[][] makePollingData() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final int ticks : Arrays.asList(1, 5, 10) ) {
            for ( final int poll : Arrays.asList(10, 100) ) {
                tests.add(new Object[]{poll, ticks});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test
    public void testPeriodUpdateNano() {
        final ProgressMeter meter = new TestingProgressMeter(10);
        final long currentTime = meter.getRuntimeInNanoseconds();
        meter.updateElapsedTimeInNanoseconds();
        Assert.assertTrue( meter.getRuntimeInNanosecondsUpdatedPeriodically() > currentTime, "Updating the periodic runtime failed" );
    }

    @Test(dataProvider = "PollingData", invocationCount = 10, successPercentage = 90, enabled = false)
    public void testProgressMeterDaemon(final long poll, final int ticks) throws InterruptedException {
        final TestingProgressMeter meter = new TestingProgressMeter(poll);
        final ProgressMeterDaemon daemon = meter.getProgressMeterDaemon();

        Assert.assertTrue(daemon.isDaemon());

        Assert.assertFalse(daemon.isDone());
        Thread.sleep(ticks * poll);
        Assert.assertFalse(daemon.isDone());

        daemon.done();
        Assert.assertTrue(daemon.isDone());

        // wait for the thread to actually finish
        daemon.join();

        Assert.assertTrue(meter.progressCalls.size() >= 1,
                "Expected at least one progress update call from daemon thread, but only got " + meter.progressCalls.size() + " with exact calls " + meter.progressCalls);

        final int tolerance = (int)Math.ceil(0.8 * meter.progressCalls.size());
        Assert.assertTrue(Math.abs(meter.progressCalls.size() - ticks) <= tolerance,
                "Expected " + ticks + " progress calls from daemon thread, but got " + meter.progressCalls.size() + " and a tolerance of only " + tolerance);

        Assert.assertTrue(meter.getRuntimeInNanosecondsUpdatedPeriodically() > 0, "Daemon should have updated our periodic runtime");
    }
}
