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

import org.apache.log4j.Logger;

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
 * Create threads, collecting statistics about their running state over time
 *
 * Uses a ThreadMXBean to capture info via ThreadInfo
 *
 * User: depristo
 * Date: 8/14/12
 * Time: 8:47 AM
 */
public class StateMonitoringThreadFactory implements ThreadFactory  {
    protected static final boolean DEBUG = false;
    private static Logger logger = Logger.getLogger(StateMonitoringThreadFactory.class);
    public static final List<Thread.State> TRACKED_STATES = Arrays.asList(Thread.State.BLOCKED, Thread.State.RUNNABLE, Thread.State.WAITING);

    final int threadsToCreate;
    final List<Thread> threads;
    final EnumMap<Thread.State, Long> times = new EnumMap<Thread.State, Long>(Thread.State.class);
    final ThreadMXBean bean;
    final CountDownLatch activeThreads;

    public StateMonitoringThreadFactory(final int threadsToCreate) {
        if ( threadsToCreate <= 0 ) throw new IllegalArgumentException("threadsToCreate <= 0: " + threadsToCreate);

        this.threadsToCreate = threadsToCreate;
        threads = new ArrayList<Thread>(threadsToCreate);
        for ( final Thread.State state : Thread.State.values() )
            times.put(state, 0l);
        bean = ManagementFactory.getThreadMXBean();
        bean.setThreadContentionMonitoringEnabled(true);
        bean.setThreadCpuTimeEnabled(true);
        activeThreads = new CountDownLatch(threadsToCreate);
    }

    public synchronized long getStateTime(final Thread.State state) {
        return times.get(state);
    }

    public synchronized long getTotalTime() {
        long total = 0;
        for ( final long time : times.values() )
            total += time;
        return total;
    }

    public synchronized double getStateFraction(final Thread.State state) {
        return getStateTime(state) / (1.0 * getTotalTime());
    }

    public int getNThreads() {
        return threads.size();
    }

    public void waitForAllThreadsToComplete() throws InterruptedException {
        activeThreads.await();
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

    @Override
    public synchronized Thread newThread(final Runnable runnable) {
        if ( threads.size() >= threadsToCreate )
            throw new IllegalStateException("Attempting to create more threads than allowed by constructor argument threadsToCreate " + threadsToCreate);

        final Thread myThread = new TrackingThread(runnable);
        threads.add(myThread);
        return myThread;
    }

    // TODO -- add polling capability

    private synchronized void updateThreadInfo(final Thread thread, final long runtime) {
        if ( DEBUG ) logger.warn("UpdateThreadInfo called");
        final ThreadInfo info = bean.getThreadInfo(thread.getId());
        if ( info != null ) {
            if ( DEBUG ) logger.warn("Updating thread total runtime " + runtime + " of which blocked " + info.getBlockedTime() + " and waiting " + info.getWaitedTime());
            incTimes(Thread.State.BLOCKED, info.getBlockedTime());
            incTimes(Thread.State.WAITING, info.getWaitedTime());
            incTimes(Thread.State.RUNNABLE, runtime - info.getWaitedTime() - info.getBlockedTime());
        }
    }

    private synchronized void incTimes(final Thread.State state, final long by) {
        times.put(state, times.get(state) + by);
    }

    private class TrackingThread extends Thread {
        private TrackingThread(Runnable runnable) {
            super(runnable);
        }

        @Override
        public void run() {
            final long startTime = System.currentTimeMillis();
            super.run();
            final long endTime = System.currentTimeMillis();
            if ( DEBUG ) logger.warn("  Countdown " + activeThreads.getCount() + " in thread " + Thread.currentThread().getName());
            updateThreadInfo(this, endTime - startTime);
            activeThreads.countDown();
            if ( DEBUG ) logger.warn("  -> Countdown " + activeThreads.getCount() + " in thread " + Thread.currentThread().getName());
        }
    }
}
