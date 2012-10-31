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

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.TimeUnit;

/**
 *
 */
public class MaxRuntimeIntegrationTest extends WalkerTest {
    private static final long STARTUP_TIME = TimeUnit.NANOSECONDS.convert(20, TimeUnit.SECONDS);

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
    @Test(enabled = true, dataProvider = "MaxRuntimeProvider", timeOut = 60 * 1000)
    public void testMaxRuntime(final MaxRuntimeTestProvider cfg) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + hg18Reference
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
}