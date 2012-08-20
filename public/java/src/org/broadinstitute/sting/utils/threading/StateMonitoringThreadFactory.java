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

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import org.apache.log4j.Logger;
import org.apache.log4j.Priority;
import org.broadinstitute.sting.utils.AutoFormattingTime;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ThreadFactory;

/**
 * Create activeThreads, collecting statistics about their running state over time
 *
 * Uses a ThreadMXBean to capture info via ThreadInfo
 *
 * User: depristo
 * Date: 8/14/12
 * Time: 8:47 AM
 */
@Invariant({
        "activeThreads.size() <= nThreadsToCreate",
        "countDownLatch.getCount() <= nThreadsToCreate",
        "nThreadsToCreated <= nThreadsToCreate"
})
public class StateMonitoringThreadFactory implements ThreadFactory  {
    protected static final boolean DEBUG = false;
    private static Logger logger = Logger.getLogger(StateMonitoringThreadFactory.class);
    public static final List<Thread.State> TRACKED_STATES = Arrays.asList(Thread.State.BLOCKED, Thread.State.RUNNABLE, Thread.State.WAITING);

    // todo -- it would be nice to not have to specify upfront the number of threads.
    // todo -- can we dynamically increment countDownLatch? It seems not...
    final int nThreadsToCreate;
    final List<Thread> activeThreads;
    final EnumMap<Thread.State, Long> times = new EnumMap<Thread.State, Long>(Thread.State.class);

    int nThreadsToCreated = 0;

    /**
     * The bean used to get the thread info about blocked and waiting times
     */
    final ThreadMXBean bean;

    /**
     * Counts down the number of active activeThreads whose runtime info hasn't been incorporated into
     * times.  Counts down from nThreadsToCreate to 0, at which point any code waiting
     * on the final times is freed to run.
     */
    final CountDownLatch countDownLatch;

    /**
     * Instead of RUNNABLE we want to print running.  This map goes from Thread.State names to human readable ones
     */
    final static EnumMap<Thread.State, String> PRETTY_NAMES = new EnumMap<Thread.State, String>(Thread.State.class);
    static {
        PRETTY_NAMES.put(Thread.State.RUNNABLE, "running");
        PRETTY_NAMES.put(Thread.State.BLOCKED,  "blocked");
        PRETTY_NAMES.put(Thread.State.WAITING,  "waiting");
    }

    /**
     * Create a new factory generating threads whose runtime and contention
     * behavior is tracked in this factory.
     *
     * @param nThreadsToCreate the number of threads we will create in the factory before it's considered complete
     *                         // TODO -- remove argument when we figure out how to implement this capability
     */
    public StateMonitoringThreadFactory(final int nThreadsToCreate) {
        if ( nThreadsToCreate <= 0 ) throw new IllegalArgumentException("nThreadsToCreate <= 0: " + nThreadsToCreate);

        this.nThreadsToCreate = nThreadsToCreate;
        activeThreads = new ArrayList<Thread>(nThreadsToCreate);

        // initialize times to 0
        for ( final Thread.State state : Thread.State.values() )
            times.put(state, 0l);

        // get the bean, and start tracking
        bean = ManagementFactory.getThreadMXBean();
        if ( bean.isThreadContentionMonitoringSupported() )
            bean.setThreadContentionMonitoringEnabled(true);
        else
            logger.warn("Thread contention monitoring not supported, we cannot track GATK multi-threaded efficiency");
            //bean.setThreadCpuTimeEnabled(true);

        countDownLatch = new CountDownLatch(nThreadsToCreate);
    }

    /**
     * Get the time spent in state across all threads created by this factory
     *
     * @param state on of the TRACKED_STATES
     * @return the time in milliseconds
     */
    @Ensures({"result >= 0", "TRACKED_STATES.contains(state)"})
    public synchronized long getStateTime(final Thread.State state) {
        return times.get(state);
    }

    /**
     * Get the total time spent in all states across all threads created by this factory
     *
     * @return the time in milliseconds
     */
    @Ensures({"result >= 0"})
    public synchronized long getTotalTime() {
        long total = 0;
        for ( final long time : times.values() )
            total += time;
        return total;
    }

    /**
     * Get the fraction of time spent in state across all threads created by this factory
     *
     * @return the fraction (0.0-1.0) of time spent in state over all state times of all threads
     */
    @Ensures({"result >= 0.0", "result <= 1.0", "TRACKED_STATES.contains(state)"})
    public synchronized double getStateFraction(final Thread.State state) {
        return getStateTime(state) / (1.0 * Math.max(getTotalTime(), 1));
    }

    /**
     * How many threads have been created by this factory so far?
     * @return
     */
    @Ensures("result >= 0")
    public int getNThreadsCreated() {
        return nThreadsToCreated;
    }

    public void waitForAllThreadsToComplete() throws InterruptedException {
        countDownLatch.await();
    }

    @Override
    public synchronized String toString() {
        final StringBuilder b = new StringBuilder();

        b.append("total ").append(getTotalTime()).append(" ");
        for ( final Thread.State state : TRACKED_STATES ) {
            b.append(state).append(" ").append(getStateTime(state)).append(" ");
        }

        return b.toString();
    }

    /**
     * Print usage information about threads from this factory to logger
     * with the INFO priority
     *
     * @param logger
     */
    public synchronized void printUsageInformation(final Logger logger) {
        printUsageInformation(logger, Priority.INFO);
    }

    /**
     * Print usage information about threads from this factory to logger
     * with the provided priority
     *
     * @param logger
     */
    public synchronized void printUsageInformation(final Logger logger, final Priority priority) {
        logger.log(priority, "Number of activeThreads used: " + getNThreadsCreated());
        logger.log(priority, "Total runtime " + new AutoFormattingTime(getTotalTime() / 1000.0));
        for ( final Thread.State state : TRACKED_STATES ) {
            logger.log(priority, String.format("  Fraction of time spent %s is %.2f (%s)",
                    prettyName(state), getStateFraction(state), new AutoFormattingTime(getStateTime(state) / 1000.0)));
        }
        logger.log(priority, String.format("Efficiency of multi-threading: %.2f%% of time spent doing productive work",
                getStateFraction(Thread.State.RUNNABLE) * 100));
    }

    private String prettyName(final Thread.State state) {
        return PRETTY_NAMES.get(state);
    }

    /**
     * Create a new thread from this factory
     *
     * @param runnable
     * @return
     */
    @Override
    @Ensures({
            "activeThreads.size() > old(activeThreads.size())",
            "activeThreads.contains(result)",
            "nThreadsToCreated == old(nThreadsToCreated) + 1"
    })
    public synchronized Thread newThread(final Runnable runnable) {
        if ( activeThreads.size() >= nThreadsToCreate)
            throw new IllegalStateException("Attempting to create more activeThreads than allowed by constructor argument nThreadsToCreate " + nThreadsToCreate);

        nThreadsToCreated++;
        final Thread myThread = new TrackingThread(runnable);
        activeThreads.add(myThread);
        return myThread;
    }

    /**
     * Update the information about completed thread that ran for runtime in milliseconds
     *
     * This method updates all of the key timing and tracking information in the factory so that
     * thread can be retired.  After this call the factory shouldn't have a pointer to the thread any longer
     *
     * @param thread
     * @param runtimeInMilliseconds
     */
    @Ensures({
            "activeThreads.size() < old(activeThreads.size())",
            "! activeThreads.contains(thread)",
            "getTotalTime() >= old(getTotalTime())",
            "countDownLatch.getCount() < old(countDownLatch.getCount())"
    })
    private synchronized void threadIsDone(final Thread thread, final long runtimeInMilliseconds) {
        if ( DEBUG ) logger.warn("  Countdown " + countDownLatch.getCount() + " in thread " + Thread.currentThread().getName());
        if ( DEBUG ) logger.warn("UpdateThreadInfo called");

        final ThreadInfo info = bean.getThreadInfo(thread.getId());
        if ( info != null ) {
            if ( DEBUG ) logger.warn("Updating thread total runtime " + runtimeInMilliseconds + " of which blocked " + info.getBlockedTime() + " and waiting " + info.getWaitedTime());
            incTimes(Thread.State.BLOCKED, info.getBlockedTime());
            incTimes(Thread.State.WAITING, info.getWaitedTime());
            incTimes(Thread.State.RUNNABLE, runtimeInMilliseconds - info.getWaitedTime() - info.getBlockedTime());
        }

        // remove the thread from the list of active activeThreads
        if ( ! activeThreads.remove(thread) )
            throw new IllegalStateException("Thread " + thread + " not in list of active activeThreads");

        // one less thread is live for those blocking on all activeThreads to be complete
        countDownLatch.countDown();
        if ( DEBUG ) logger.warn("  -> Countdown " + countDownLatch.getCount() + " in thread " + Thread.currentThread().getName());
    }

    /**
     * Helper function that increments the times counter by by for state
     *
     * @param state
     * @param by
     */
    private synchronized void incTimes(final Thread.State state, final long by) {
        times.put(state, times.get(state) + by);
    }

    /**
     * A wrapper around Thread that tracks the runtime of the thread and calls threadIsDone() when complete
     */
    private class TrackingThread extends Thread {
        private TrackingThread(Runnable runnable) {
            super(runnable);
        }

        @Override
        public void run() {
            final long startTime = System.currentTimeMillis();
            super.run();
            final long endTime = System.currentTimeMillis();
            threadIsDone(this, endTime - startTime);
        }
    }
}
