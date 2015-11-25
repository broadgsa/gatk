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

package org.broadinstitute.gatk.engine;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.SimpleTimer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 *
 */
public class MaxRuntimeIntegrationTest extends WalkerTest {
    public static class SleepingWalker extends LocusWalker<Integer, Integer> {
        @Output PrintStream out;

        @Argument(fullName="sleepTime",shortName="sleepTime",doc="x", required=false)
        public int sleepTime = 100;

        @Override
        public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            try {Thread.sleep(sleepTime);} catch (InterruptedException e) {};
            return 1;
        }

        @Override public Integer reduceInit() { return 0; }
        @Override public Integer reduce(Integer value, Integer sum) { return sum + value; }

        @Override
        public void onTraversalDone(Integer result) {
            out.println(result);
        }
    }

    private static final long STARTUP_TIME = TimeUnit.NANOSECONDS.convert(60, TimeUnit.SECONDS);

    private class MaxRuntimeTestProvider extends TestDataProvider {
        final long maxRuntime;
        final TimeUnit unit;

        public MaxRuntimeTestProvider(final long maxRuntime, final TimeUnit unit) {
            super(MaxRuntimeTestProvider.class);
            this.maxRuntime = maxRuntime;
            this.unit = unit;
            setName(String.format("Max runtime test : %d of %s", maxRuntime, unit));
        }

        public long expectedMaxRuntimeNano() {
            return TimeUnit.NANOSECONDS.convert(maxRuntime, unit) + STARTUP_TIME;
        }
    }

    @DataProvider(name = "MaxRuntimeProvider")
    public Object[][] makeMaxRuntimeProvider() {
        for ( final TimeUnit requestedUnits : Arrays.asList(TimeUnit.NANOSECONDS, TimeUnit.MILLISECONDS, TimeUnit.SECONDS, TimeUnit.MINUTES) )
            new MaxRuntimeTestProvider(requestedUnits.convert(30, TimeUnit.SECONDS), requestedUnits);

        return MaxRuntimeTestProvider.getTests(MaxRuntimeTestProvider.class);
    }

    //
    // Loop over errors to throw, make sure they are the errors we get back from the engine, regardless of NT type
    //
    @Test(enabled = true, dataProvider = "MaxRuntimeProvider", timeOut = 120 * 1000)
    public void testMaxRuntime(final MaxRuntimeTestProvider cfg) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestPrintReadsWalker -R " + hg18Reference
                        + " -I " + validationDataLocation + "NA12878.WEx.downsampled20x.bam -o /dev/null"
                        + " -maxRuntime " + cfg.maxRuntime + " -maxRuntimeUnits " + cfg.unit, 0,
                Collections.<String>emptyList());
        final SimpleTimer timer = new SimpleTimer().start();
        executeTest("Max runtime " + cfg, spec);
        final long actualRuntimeNano = timer.getElapsedTimeNano();

        Assert.assertTrue(actualRuntimeNano < cfg.expectedMaxRuntimeNano(),
                "Actual runtime " + TimeUnit.SECONDS.convert(actualRuntimeNano, TimeUnit.NANOSECONDS)
                        + " exceeded max. tolerated runtime " + TimeUnit.SECONDS.convert(cfg.expectedMaxRuntimeNano(), TimeUnit.NANOSECONDS)
                        + " given requested runtime " + cfg.maxRuntime + " " + cfg.unit);
    }

    @DataProvider(name = "SubshardProvider")
    public Object[][] makeSubshardProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{10});
        tests.add(new Object[]{100});
        tests.add(new Object[]{500});
        tests.add(new Object[]{1000});
        tests.add(new Object[]{2000});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "SubshardProvider", timeOut = 120 * 1000)
    public void testSubshardTimeout(final int sleepTime) throws Exception {
        final int maxRuntime = 5000;

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T SleepingWalker -R " + b37KGReference
                        + " -I " + privateTestDir + "NA12878.100kb.BQSRv2.example.bam -o %s"
                        + " -maxRuntime " + maxRuntime + " -maxRuntimeUnits MILLISECONDS -sleepTime " + sleepTime, 1,
                Collections.singletonList(""));
        final File result = executeTest("Subshard max runtime ", spec).getFirst().get(0);
        final int cycle = Integer.valueOf(new BufferedReader(new FileReader(result)).readLine());

        final int maxCycles = (int)Math.ceil((maxRuntime * 5) / sleepTime);
        logger.warn(String.format("Max cycles %d saw %d in file %s with sleepTime %d and maxRuntime %d", maxCycles, cycle, result, sleepTime, maxRuntime));
        Assert.assertTrue(cycle < maxCycles, "Too many cycles seen -- saw " + cycle + " in file " + result + " but max should have been " + maxCycles);
    }
}