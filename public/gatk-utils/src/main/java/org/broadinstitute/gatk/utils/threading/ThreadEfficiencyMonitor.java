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

package org.broadinstitute.gatk.utils.threading;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.apache.log4j.Priority;
import org.broadinstitute.gatk.utils.AutoFormattingTime;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.util.EnumMap;
import java.util.concurrent.TimeUnit;

/**
 * Uses an MXBean to monitor thread efficiency
 *
 * Once the monitor is created, calls to threadIsDone() can be used to add information
 * about the efficiency of the provided thread to this monitor.
 *
 * Provides simple print() for displaying efficiency information to a logger
 *
 * User: depristo
 * Date: 8/22/12
 * Time: 10:48 AM
 */
@Invariant({"nThreadsAnalyzed >= 0"})
public class ThreadEfficiencyMonitor {
    protected static final boolean DEBUG = false;
    protected static Logger logger = Logger.getLogger(EfficiencyMonitoringThreadFactory.class);
    final EnumMap<State, Long> times = new EnumMap<State, Long>(State.class);

    /**
     * The number of threads we've included in our efficiency monitoring
     */
    int nThreadsAnalyzed = 0;

    /**
     * The bean used to get the thread info about blocked and waiting times
     */
    final ThreadMXBean bean;

    public ThreadEfficiencyMonitor() {
        bean = ManagementFactory.getThreadMXBean();

        // get the bean, and start tracking
        if ( bean.isThreadContentionMonitoringSupported() )
            bean.setThreadContentionMonitoringEnabled(true);
        else
            logger.warn("Thread contention monitoring not supported, we cannot track GATK multi-threaded efficiency");
        //bean.setThreadCpuTimeEnabled(true);

        if ( bean.isThreadCpuTimeSupported() )
            bean.setThreadCpuTimeEnabled(true);
        else
            logger.warn("Thread CPU monitoring not supported, we cannot track GATK multi-threaded efficiency");

        // initialize times to 0
        for ( final State state : State.values() )
            times.put(state, 0l);
    }

    private static long nanoToMilli(final long timeInNano) {
        return TimeUnit.NANOSECONDS.toMillis(timeInNano);
    }

    /**
     * Get the time spent in state across all threads created by this factory
     *
     * @param state to get information about
     * @return the time in milliseconds
     */
    @Ensures({"result >= 0"})
    public synchronized long getStateTime(final State state) {
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
     * @return the percentage (0.0-100.0) of time spent in state over all state times of all threads
     */
    @Ensures({"result >= 0.0", "result <= 100.0"})
    public synchronized double getStatePercent(final State state) {
        return (100.0 * getStateTime(state)) / Math.max(getTotalTime(), 1);
    }

    public int getnThreadsAnalyzed() {
        return nThreadsAnalyzed;
    }

    @Override
    public synchronized String toString() {
        final StringBuilder b = new StringBuilder();

        b.append("total ").append(getTotalTime()).append(" ");
        for ( final State state : State.values() ) {
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
        logger.debug("Number of threads monitored: " + getnThreadsAnalyzed());
        logger.debug("Total runtime " + new AutoFormattingTime(TimeUnit.MILLISECONDS.toNanos(getTotalTime())));
        for ( final State state : State.values() ) {
            logger.debug(String.format("\tPercent of time spent %s is %.2f", state.getUserFriendlyName(), getStatePercent(state)));
        }
        logger.log(priority, String.format("CPU      efficiency : %6.2f%% of time spent %s", getStatePercent(State.USER_CPU), State.USER_CPU.getUserFriendlyName()));
        logger.log(priority, String.format("Walker inefficiency : %6.2f%% of time spent %s", getStatePercent(State.BLOCKING), State.BLOCKING.getUserFriendlyName()));
        logger.log(priority, String.format("I/O    inefficiency : %6.2f%% of time spent %s", getStatePercent(State.WAITING_FOR_IO), State.WAITING_FOR_IO.getUserFriendlyName()));
        logger.log(priority, String.format("Thread inefficiency : %6.2f%% of time spent %s", getStatePercent(State.WAITING), State.WAITING.getUserFriendlyName()));
    }

    /**
     * Update the information about completed thread that ran for runtime in milliseconds
     *
     * This method updates all of the key timing and tracking information in the factory so that
     * thread can be retired.  After this call the factory shouldn't have a pointer to the thread any longer
     *
     * @param thread the thread whose information we are updating
     */
    @Ensures({
            "getTotalTime() >= old(getTotalTime())"
    })
    public synchronized void threadIsDone(final Thread thread) {
        nThreadsAnalyzed++;

        if ( DEBUG ) logger.warn("UpdateThreadInfo called");

        final long threadID = thread.getId();
        final ThreadInfo info = bean.getThreadInfo(thread.getId());
        final long totalTimeNano = bean.getThreadCpuTime(threadID);
        final long userTimeNano = bean.getThreadUserTime(threadID);
        final long systemTimeNano = totalTimeNano - userTimeNano;
        final long userTimeInMilliseconds = nanoToMilli(userTimeNano);
        final long systemTimeInMilliseconds = nanoToMilli(systemTimeNano);

        if ( info != null ) {
            if ( DEBUG ) logger.warn("Updating thread with user runtime " + userTimeInMilliseconds + " and system runtime " + systemTimeInMilliseconds + " of which blocked " + info.getBlockedTime() + " and waiting " + info.getWaitedTime());
            incTimes(State.BLOCKING, info.getBlockedTime());
            incTimes(State.WAITING, info.getWaitedTime());
            incTimes(State.USER_CPU, userTimeInMilliseconds);
            incTimes(State.WAITING_FOR_IO, systemTimeInMilliseconds);
        }
    }

    /**
     * Helper function that increments the times counter by by for state
     *
     * @param state
     * @param by
     */
    @Requires({"state != null", "by >= 0"})
    @Ensures("getTotalTime() == old(getTotalTime()) + by")
    private synchronized void incTimes(final State state, final long by) {
        times.put(state, times.get(state) + by);
    }

    public enum State {
        BLOCKING("blocking on synchronized data structures"),
        WAITING("waiting on some other thread"),
        USER_CPU("doing productive CPU work"),
        WAITING_FOR_IO("waiting for I/O");

        private final String userFriendlyName;

        private State(String userFriendlyName) {
            this.userFriendlyName = userFriendlyName;
        }

        public String getUserFriendlyName() {
            return userFriendlyName;
        }
    }
}
