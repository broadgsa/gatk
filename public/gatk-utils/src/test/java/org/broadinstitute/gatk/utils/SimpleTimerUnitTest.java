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

package org.broadinstitute.gatk.utils;

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.lang.reflect.Field;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class SimpleTimerUnitTest extends BaseTest {
    private final static String NAME = "unit.test.timer";

    @Test
    public void testSimpleTimer() {
        SimpleTimer t = new SimpleTimer(NAME);
        Assert.assertEquals(t.getName(), NAME, "Name is not the provided one");
        Assert.assertFalse(t.isRunning(), "Initial state of the timer is running");
        Assert.assertEquals(t.getElapsedTime(), 0.0, "New timer elapsed time should be 0");
        Assert.assertEquals(t.getElapsedTimeNano(), 0l, "New timer elapsed time nano should be 0");

        t.start();
        Assert.assertTrue(t.isRunning(), "Started timer isn't running");
        Assert.assertTrue(t.getElapsedTime() >= 0.0, "Elapsed time should be >= 0");
        Assert.assertTrue(t.getElapsedTimeNano() >= 0.0, "Elapsed time nano should be >= 0");
        long n1 = t.getElapsedTimeNano();
        double t1 = t.getElapsedTime();
        idleLoop(); // idle loop to wait a tiny bit of time
        long n2 = t.getElapsedTimeNano();
        double t2 = t.getElapsedTime();
        Assert.assertTrue(t2 >= t1, "T2 >= T1 for a running time");
        Assert.assertTrue(n2 >= n1, "T2 >= T1 nano for a running time");

        t.stop();
        Assert.assertFalse(t.isRunning(), "Stopped timer still running");
        long n3 = t.getElapsedTimeNano();
        double t3 = t.getElapsedTime();
        idleLoop(); // idle loop to wait a tiny bit of time
        double t4 = t.getElapsedTime();
        long n4 = t.getElapsedTimeNano();
        Assert.assertTrue(t4 == t3, "Elapsed times for two calls of stop timer not the same");
        Assert.assertTrue(n4 == n3, "Elapsed times for two calls of stop timer not the same");

        t.restart();
        idleLoop(); // idle loop to wait a tiny bit of time
        double t5 = t.getElapsedTime();
        long n5 = t.getElapsedTimeNano();
        Assert.assertTrue(t.isRunning(), "Restarted timer should be running");
        idleLoop(); // idle loop to wait a tiny bit of time
        double t6 = t.getElapsedTime();
        long n6 = t.getElapsedTimeNano();
        Assert.assertTrue(t5 >= t4, "Restarted timer elapsed time should be after elapsed time preceding the restart");
        Assert.assertTrue(t6 >= t5, "Second elapsed time not after the first in restarted timer");
        Assert.assertTrue(n5 >= n4, "Restarted timer elapsed time nano should be after elapsed time preceding the restart");
        Assert.assertTrue(n6 >= n5, "Second elapsed time nano not after the first in restarted timer");

        final List<Double> secondTimes = Arrays.asList(t1, t2, t3, t4, t5, t6);
        final List<Long> nanoTimes     = Arrays.asList(n1, n2, n3, n4, n5, n6);
        for ( int i = 0; i < nanoTimes.size(); i++ )
            Assert.assertEquals(
                    SimpleTimer.nanoToSecondsAsDouble(nanoTimes.get(i)),
                    secondTimes.get(i), 1e-1, "Nanosecond and second timer disagree");
    }

    @Test
    public void testNanoResolution() {
        SimpleTimer t = new SimpleTimer(NAME);

        // test the nanosecond resolution
        long n7 = t.currentTimeNano();
        int sum = 0;
        for ( int i = 0; i < 100; i++) sum += i;
        long n8 = t.currentTimeNano();
        final long delta = n8 - n7;
        final long oneMilliInNano = TimeUnit.MILLISECONDS.toNanos(1);
        logger.warn("nanoTime before nano operation " + n7);
        logger.warn("nanoTime after nano operation of summing 100 ints " + n8 + ", sum = " + sum + " time delta " + delta + " vs. 1 millsecond in nano " + oneMilliInNano);
        Assert.assertTrue(n8 > n7, "SimpleTimer doesn't appear to have nanoSecond resolution: n8 " + n8 + " <= n7 " + n7);
        Assert.assertTrue(delta < oneMilliInNano,
                "SimpleTimer doesn't appear to have nanoSecond resolution: time delta is " + delta + " vs 1 millisecond in nano " + oneMilliInNano);
    }

    @Test
    public void testMeaningfulTimes() {
        SimpleTimer t = new SimpleTimer(NAME);

        t.start();
        for ( int i = 0; i < 100; i++ ) ;
        long nano = t.getElapsedTimeNano();
        double secs = t.getElapsedTime();

        Assert.assertTrue(secs > 0, "Seconds timer doesn't appear to count properly: elapsed time is " + secs);
        Assert.assertTrue(secs < 0.01, "Fast operation said to take longer than 10 milliseconds: elapsed time in seconds " + secs);

        Assert.assertTrue(nano > 0, "Nanosecond timer doesn't appear to count properly: elapsed time is " + nano);
        final long maxTimeInMicro = 10000;
        final long maxTimeInNano = TimeUnit.MICROSECONDS.toNanos(maxTimeInMicro);
        Assert.assertTrue(nano < maxTimeInNano, "Fast operation said to take longer than " + maxTimeInMicro + " microseconds: elapsed time in nano " + nano + " micro " + TimeUnit.NANOSECONDS.toMicros(nano));
    }

    @Test
    public void testCheckpointRestart() throws Exception {
        SimpleTimer t = new SimpleTimer(NAME);
        
        final Field offsetField = t.getClass().getDeclaredField("nanoTimeOffset");
        offsetField.setAccessible(true);
        long offset = ((Long) offsetField.get(t)).longValue();

        t.start();
        idleLoop();
        // Make it as if clock has jumped into the past
        offsetField.set(t, offset + TimeUnit.SECONDS.toNanos(10));
        t.stop();
        offset = ((Long) offsetField.get(t)).longValue();
        Assert.assertEquals(t.getElapsedTime(), 0.0, "Time over restart is not zero.");

        t.start();
        idleLoop();
        t.stop();
        offset = ((Long) offsetField.get(t)).longValue();
        double elapsed = t.getElapsedTime();
        Assert.assertTrue(elapsed >= 0.0, "Elapsed time is zero.");
        t.restart();
        // Make the clock jump again by just a little
        offsetField.set(t, offset + TimeUnit.SECONDS.toNanos(1));
        idleLoop();
        t.stop();
        offset = ((Long) offsetField.get(t)).longValue();
        Assert.assertTrue(t.getElapsedTime() > elapsed, "Small clock drift causing reset.");
        elapsed = t.getElapsedTime();
        // Now a bigger jump, into the future this time.
        t.restart();
        // Make the clock jump again by a lot
        offsetField.set(t, offset - TimeUnit.SECONDS.toNanos(10));
        t.stop();
        Assert.assertEquals(t.getElapsedTime(), elapsed, "Time added over checkpoint/restart.");

        // Test without stopping
        t.start();
        offset = ((Long) offsetField.get(t)).longValue();
        // Make it as if clock has jumped into the past
        offsetField.set(t, offset + TimeUnit.SECONDS.toNanos(10));       
        Assert.assertEquals(t.getElapsedTime(), 0.0, "Elapsed time after C/R is not zero.");
        idleLoop();
        Assert.assertTrue(t.getElapsedTime() > 0.0, "Elapsed time zero after re-sync.");

    }

    private static void idleLoop() {
        for ( int i = 0; i < 100000; i++ ) ; // idle loop to wait a tiny bit of time
    }
}