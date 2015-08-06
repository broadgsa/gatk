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

package org.broadinstitute.gatk.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.MultiThreadedErrorTracker;

import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Reducer supporting multi-threaded reduce of the map/reduce.
 *
 * reduceAsMuchAsPossible is the key function.  Multiple threads can call into this, providing
 * the map results queue, and this class accumulates the result of calling reduce
 * on the maps objects.  reduceAsMuchAsPossible isn't directly synchronized, but manages multi-threading
 * directly with a lock.  Threads can request either to block on the reduce call until it can be
 * executed, or immediately exit if the lock isn't available.  That allows multi-threaded users
 * to avoid piling up waiting to reduce while one thread is reducing.  They can instead immediately
 * leave to go do something else productive
 *
 * @author depristo
 * @since 2012
 */
class Reducer<MapType, ReduceType> {
    private final static Logger logger = Logger.getLogger(Reducer.class);

    /**
     * The reduce function to execute
     */
    private final NSReduceFunction<MapType, ReduceType> reduce;

    /**
     * Used to communicate errors to the outer master thread
     */
    private final MultiThreadedErrorTracker errorTracker;

    /**
     * Lock used to protect the call reduceAsMuchAsPossible from race conditions
     */
    private final Lock reduceLock = new ReentrantLock();

    /**
     * The sum of the reduce function applied to all MapResults.  After this Reducer
     * is done sum contains the final reduce result.
     */
    ReduceType sum;

    /**
     * Create a new Reducer that will apply the reduce function with initialSum value
     * to values via reduceAsMuchAsPossible, timing the reduce function call costs with
     * reduceTimer
     *
     * @param reduce the reduce function to apply
     * @param initialSum the initial reduce sum
     */
    public Reducer(final NSReduceFunction<MapType, ReduceType> reduce,
                   final MultiThreadedErrorTracker errorTracker,
                   final ReduceType initialSum) {
        if ( errorTracker == null ) throw new IllegalArgumentException("Error tracker cannot be null");
        if ( reduce == null ) throw new IllegalArgumentException("Reduce function cannot be null");

        this.errorTracker = errorTracker;
        this.reduce = reduce;
        this.sum = initialSum;
    }

    /**
     * Reduce as much data as possible in mapResultQueue, returning the number of reduce calls completed
     *
     * As much as possible is defined as all of the MapResults in the queue are in order starting from the
     * numSubmittedJobs we reduced previously, up to the either the queue being empty or where the next MapResult
     * doesn't have JobID == prevJobID + 1.
     *
     * @param mapResultQueue a queue of MapResults in jobID order
     * @return the number of reduces run, from 0 >
     * @throws InterruptedException
     */
    @Ensures("result >= 0")
    public int reduceAsMuchAsPossible(final MapResultsQueue<MapType> mapResultQueue, final boolean waitForLock) {
        if ( mapResultQueue == null ) throw new IllegalArgumentException("mapResultQueue cannot be null");
        int nReducesNow = 0;

        final boolean haveLock = acquireReduceLock(waitForLock);
        try {
            if ( haveLock ) {
                while ( mapResultQueue.nextValueIsAvailable() ) {
                    final MapResult<MapType> result = mapResultQueue.take();

                    if ( ! result.isEOFMarker() ) {
                        nReducesNow++;

                        // apply reduce, keeping track of sum
                        sum = reduce.apply(result.getValue(), sum);
                    }
                }
            }
        } catch (Exception ex) {
            errorTracker.notifyOfError(ex);
        } finally {
            if ( haveLock ) // if we acquired the lock, unlock it
                releaseReduceLock();
        }

        return nReducesNow;
    }

    /**
     * Acquire the reduce lock, either returning immediately if not possible or blocking until the lock is available
     *
     * @param blockUntilAvailable if true, we will block until the lock is available, otherwise we return immediately
     *                            without acquiring the lock
     * @return true if the lock has been acquired, false otherwise
     */
    protected boolean acquireReduceLock(final boolean blockUntilAvailable) {
        if ( blockUntilAvailable ) {
            reduceLock.lock();
            return true;
        } else {
            return reduceLock.tryLock();
        }
    }

    /**
     * Free the reduce lock.
     *
     * Assumes that the invoking thread actually previously acquired the lock (it's a problem if not).
     */
    protected void releaseReduceLock() {
        reduceLock.unlock();
    }

    /**
     * Get the current reduce result resulting from applying reduce(...) to all MapResult elements.
     *
     * Note that this method cannot know if future reduce calls are coming in.  So it simply gets
     * the current reduce result.  It is up to the caller to know whether the returned value is
     * a partial result, or the full final value
     *
     * @return the total reduce result across all jobs
     */
    public ReduceType getReduceResult() {
        return sum;
    }
}
