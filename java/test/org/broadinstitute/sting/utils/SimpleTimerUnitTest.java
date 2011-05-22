package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;

public class SimpleTimerUnitTest extends BaseTest {
    private final static String NAME = "unit.test.timer";

    @Test
    public void testSimpleTimer() {
        SimpleTimer t = new SimpleTimer(NAME);
        Assert.assertEquals(t.getName(), NAME, "Name is not the provided one");
        Assert.assertFalse(t.isRunning(), "Initial state of the timer is running");
        Assert.assertEquals(t.getElapsedTime(), 0.0, "New timer elapsed time should be 0");

        t.start();
        Assert.assertTrue(t.isRunning(), "Started timer isn't running");
        Assert.assertTrue(t.getElapsedTime() >= 0.0, "Elapsed time should be >= 0");
        double t1 = t.getElapsedTime();
        idleLoop(); // idle loop to wait a tiny bit of time
        double t2 = t.getElapsedTime();
        Assert.assertTrue(t2 >= t1, "T2 >= T1 for a running time");

        t.stop();
        Assert.assertFalse(t.isRunning(), "Stopped timer still running");
        double t3 = t.getElapsedTime();
        idleLoop(); // idle loop to wait a tiny bit of time
        double t4 = t.getElapsedTime();
        Assert.assertTrue(t4 == t3, "Elapsed times for two calls of stop timer not the same");

        t.restart();
        idleLoop(); // idle loop to wait a tiny bit of time
        double t5 = t.getElapsedTime();
        Assert.assertTrue(t.isRunning(), "Restarted timer should be running");
        idleLoop(); // idle loop to wait a tiny bit of time
        double t6 = t.getElapsedTime();
        Assert.assertTrue(t5 >= t4, "Restarted timer elapsed time should be after elapsed time preceding the restart");
        Assert.assertTrue(t6 >= t5, "Second elapsed time not after the first in restarted timer");

        t.stop().start();
        Assert.assertTrue(t.isRunning(), "second started timer isn't running");
        Assert.assertTrue(t.getElapsedTime() >= 0.0, "elapsed time should have been reset");
        Assert.assertTrue(t.getElapsedTime() < t6, "elapsed time isn't less than time before start call"); // we should have effective no elapsed time
    }

    private final static void idleLoop() {
        for ( int i = 0; i < 100000; i++ ) ; // idle loop to wait a tiny bit of time
    }
}