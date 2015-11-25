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
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;

/**
 * Creates threads that automatically monitor their efficiency via the parent ThreadEfficiencyMonitor
 *
 * User: depristo
 * Date: 8/14/12
 * Time: 8:47 AM
 */
@Invariant({
        "activeThreads.size() <= nThreadsToCreate",
        "countDownLatch.getCount() <= nThreadsToCreate",
        "nThreadsCreated <= nThreadsToCreate"
})
public class EfficiencyMonitoringThreadFactory extends ThreadEfficiencyMonitor implements ThreadFactory  {
    final int nThreadsToCreate;
    final List<Thread> activeThreads;

    int nThreadsCreated = 0;

    /**
     * Counts down the number of active activeThreads whose runtime info hasn't been incorporated into
     * times.  Counts down from nThreadsToCreate to 0, at which point any code waiting
     * on the final times is freed to run.
     */
    final CountDownLatch countDownLatch;

    /**
     * Create a new factory generating threads whose runtime and contention
     * behavior is tracked in this factory.
     *
     * @param nThreadsToCreate the number of threads we will create in the factory before it's considered complete
     */
    public EfficiencyMonitoringThreadFactory(final int nThreadsToCreate) {
        super();
        if ( nThreadsToCreate <= 0 ) throw new IllegalArgumentException("nThreadsToCreate <= 0: " + nThreadsToCreate);

        this.nThreadsToCreate = nThreadsToCreate;
        activeThreads = new ArrayList<Thread>(nThreadsToCreate);
        countDownLatch = new CountDownLatch(nThreadsToCreate);
    }

    /**
     * How many threads have been created by this factory so far?
     * @return
     */
    @Ensures("result >= 0")
    public int getNThreadsCreated() {
        return nThreadsCreated;
    }

    /**
     * Only useful for testing, so that we can wait for all of the threads in the factory to complete running
     *
     * @throws InterruptedException
     */
    protected void waitForAllThreadsToComplete() throws InterruptedException {
        countDownLatch.await();
    }

    @Ensures({
            "activeThreads.size() <= old(activeThreads.size())",
            "! activeThreads.contains(thread)",
            "countDownLatch.getCount() <= old(countDownLatch.getCount())"
    })
    @Override
    public synchronized void threadIsDone(final Thread thread) {
        nThreadsAnalyzed++;

        if ( DEBUG ) logger.warn("  Countdown " + countDownLatch.getCount() + " in thread " + Thread.currentThread().getName());

        super.threadIsDone(thread);

        // remove the thread from the list of active activeThreads, if it's in there, and decrement the countdown latch
        if ( activeThreads.remove(thread) ) {
            // one less thread is live for those blocking on all activeThreads to be complete
            countDownLatch.countDown();
            if ( DEBUG ) logger.warn("  -> Countdown " + countDownLatch.getCount() + " in thread " + Thread.currentThread().getName());
        }
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
            "nThreadsCreated == old(nThreadsCreated) + 1"
    })
    public synchronized Thread newThread(final Runnable runnable) {
        if ( activeThreads.size() >= nThreadsToCreate)
            throw new IllegalStateException("Attempting to create more activeThreads than allowed by constructor argument nThreadsToCreate " + nThreadsToCreate);

        nThreadsCreated++;
        final Thread myThread = new TrackingThread(runnable);
        activeThreads.add(myThread);
        return myThread;
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
            super.run();
            threadIsDone(this);
        }
    }
}
