/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Priority;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;

/**
 * Tests for the state monitoring thread factory.
 */
public class StateMonitoringThreadFactoryUnitTest extends BaseTest {
    // the duration of the tests -- 100 ms is tolerable given the number of tests we are doing
    private final static long THREAD_TARGET_DURATION_IN_MILLISECOND = 100;
    final static Object GLOBAL_LOCK = new Object();

    private class StateTest extends TestDataProvider {
        private final double TOLERANCE = 0.1; // willing to tolerate a 10% error

        final List<Thread.State> statesForThreads;

        public StateTest(final List<Thread.State> statesForThreads) {
            super(StateTest.class);
            this.statesForThreads = statesForThreads;
            setName("StateTest " + Utils.join(",", statesForThreads));
        }

        public List<Thread.State> getStatesForThreads() {
            return statesForThreads;
        }

        public int getNStates() { return statesForThreads.size(); }

        public double maxStateFraction(final Thread.State state) { return fraction(state) + TOLERANCE; }
        public double minStateFraction(final Thread.State state) { return fraction(state) - TOLERANCE; }

        private double fraction(final Thread.State state) {
            return Collections.frequency(statesForThreads, state) / (1.0 * statesForThreads.size());
        }
    }

    /**
     * Test helper threading class that puts the thread into RUNNING, BLOCKED, or WAITING state as
     * requested for input argument
     */
    private static class StateTestThread implements Callable<Double> {
        private final Thread.State stateToImplement;

        private StateTestThread(final Thread.State stateToImplement) {
            if ( ! StateMonitoringThreadFactory.TRACKED_STATES.contains(stateToImplement) )
                throw new IllegalArgumentException("Unexpected state " + stateToImplement);
            this.stateToImplement = stateToImplement;
        }

        @Override
        public Double call() throws Exception {
            switch ( stateToImplement ) {
                case RUNNABLE:
                    // do some work until we get to THREAD_TARGET_DURATION_IN_MILLISECOND
                    double sum = 0.0;
                    final long startTime = System.currentTimeMillis();
                    for ( int i = 1; System.currentTimeMillis() - startTime < (THREAD_TARGET_DURATION_IN_MILLISECOND - 1); i++ ) {
                        sum += Math.log10(i);
                    }
                    return sum;
                case WAITING:
                    Thread.currentThread().sleep(THREAD_TARGET_DURATION_IN_MILLISECOND);
                    return 0.0;
                case BLOCKED:
                    if ( StateMonitoringThreadFactory.DEBUG ) logger.warn("Blocking...");
                    synchronized (GLOBAL_LOCK) {
                        // the GLOBAL_LOCK must be held by the unit test itself for this to properly block
                        if ( StateMonitoringThreadFactory.DEBUG ) logger.warn("  ... done blocking");
                    }
                    return 0.0;
                default:
                    throw new ReviewedStingException("Unexpected thread test state " + stateToImplement);
            }
        }
    }

    @DataProvider(name = "StateTest")
    public Object[][] createStateTest() {
        for ( final int nThreads : Arrays.asList(1, 2, 3, 4) ) {
            for (final List<Thread.State> states : Utils.makeCombinations(StateMonitoringThreadFactory.TRACKED_STATES, nThreads) ) {
                //if ( Collections.frequency(states, Thread.State.BLOCKED) > 0)
                    new StateTest(states);
            }
        }

        return StateTest.getTests(StateTest.class);
    }

    @Test(enabled = true, dataProvider = "StateTest")
    public void testStateTest(final StateTest test) throws InterruptedException {
        // allows us to test blocking
        final StateMonitoringThreadFactory factory = new StateMonitoringThreadFactory(test.getNStates());
        final ExecutorService threadPool = Executors.newFixedThreadPool(test.getNStates(), factory);

        logger.warn("Running " + test);
        synchronized (GLOBAL_LOCK) {
            //logger.warn("  Have lock");
            for ( final Thread.State threadToRunState : test.getStatesForThreads() )
            threadPool.submit(new StateTestThread(threadToRunState));

            // lock has to be here for the whole running of the activeThreads but end before the sleep so the blocked activeThreads
            // can block for their allotted time
            threadPool.shutdown();
            Thread.sleep(THREAD_TARGET_DURATION_IN_MILLISECOND);
        }
        //logger.warn("  Releasing lock");
        threadPool.awaitTermination(10, TimeUnit.SECONDS);
        //logger.warn("  done awaiting termination");
        //logger.warn("  waiting for all activeThreads to complete");
        factory.waitForAllThreadsToComplete();
        //logger.warn("  done waiting for activeThreads");

        // make sure we counted everything properly
        final long totalTime = factory.getTotalTime();
        final long minTime = (THREAD_TARGET_DURATION_IN_MILLISECOND - 10) * test.getNStates();
        final long maxTime = (THREAD_TARGET_DURATION_IN_MILLISECOND + 10) * test.getNStates();
        //logger.warn("Testing total time");
        Assert.assertTrue(totalTime >= minTime, "Factory results not properly accumulated: totalTime = " + totalTime + " < minTime = " + minTime);
        Assert.assertTrue(totalTime <= maxTime, "Factory results not properly accumulated: totalTime = " + totalTime + " > maxTime = " + maxTime);

        for (final Thread.State state : StateMonitoringThreadFactory.TRACKED_STATES ) {
            final double min = test.minStateFraction(state);
            final double max = test.maxStateFraction(state);
            final double obs = factory.getStateFraction(state);
//            logger.warn("  Checking " + state
//                    + " min " + String.format("%.2f", min)
//                    + " max " + String.format("%.2f", max)
//                    + " obs " + String.format("%.2f", obs)
//                    + " factor = " + factory);
            Assert.assertTrue(obs >= min, "Too little time spent in state " + state + " obs " + obs + " min " + min);
            Assert.assertTrue(obs <= max, "Too much time spent in state " + state + " obs " + obs + " max " + min);
        }

        // we actually ran the expected number of activeThreads
        Assert.assertEquals(factory.getNThreadsCreated(), test.getNStates());

        // should be called to ensure we don't format / NPE on output
        factory.printUsageInformation(logger, Priority.INFO);
    }
}