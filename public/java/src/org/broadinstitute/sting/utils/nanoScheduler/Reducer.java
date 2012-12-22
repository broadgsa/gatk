package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MultiThreadedErrorTracker;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Reducer supporting two-threaded reduce of the map/reduce.
 *
 * The first thread, using the reduceAsMuchAsPossible function, actually reduces the data
 * as it arrives in the blockingQueue.
 *
 * The second thread, using the waitForFinalReduce, can block on this data structure
 * until that all jobs have arrived and been reduced.
 *
 * The key function for communication here is setTotalJobCount(), which the thread that submits
 * jobs that enqueue MapResults into the blocking queue must call ONCE to tell the
 * Reducer the total number of jobs that have been submitted for map.  When numOfSubmittedJobs
 * have been processed, this class frees a latch that allows thread blocked on waitForFinalReduce to proceed.
 *
 * This thread reads from mapResultsQueue until the poison EOF object arrives.  At each
 * stage is calls reduce(value, sum).  The blocking mapResultQueue ensures that the
 * queue waits until the mapResultQueue has a value to take. Then, it gets and waits
 * until the map result Future has a value.
 */
class Reducer<MapType, ReduceType> {
    private final static Logger logger = Logger.getLogger(Reducer.class);
    private final static int UNSET_NUM_SUBMITTED_JOBS = -2;

    private final CountDownLatch countDownLatch = new CountDownLatch(1);
    private final NSReduceFunction<MapType, ReduceType> reduce;
    private final MultiThreadedErrorTracker errorTracker;
    private final Lock reduceLock = new ReentrantLock();

    /**
     * The sum of the reduce function applied to all MapResults.  After this Reducer
     * is done sum contains the final reduce result.
     */
    ReduceType sum;

    int numSubmittedJobs = UNSET_NUM_SUBMITTED_JOBS; // not yet set

    /**
     * A counter keeping track of the number of jobs we're reduced
     */
    int numJobsReduced = 0;

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

                    numJobsReduced++;
                    maybeReleaseLatch();
                }
            }
        } catch (Exception ex) {
            errorTracker.notifyOfError(ex);
            countDownLatch.countDown();
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
     * release the latch if appropriate
     *
     * Appropriate means we've seen the last job, or there's only a single job id
     */
    private void maybeReleaseLatch() {
        if ( numJobsReduced == numSubmittedJobs ) {
            // either we've already seen the last one prevJobID == numSubmittedJobs or
            // the last job ID is -1, meaning that no jobs were ever submitted
            countDownLatch.countDown();
        }
    }

    /**
     * For testing only
     *
     * @return true if latch is released
     */
    protected synchronized boolean latchIsReleased() {
        return countDownLatch.getCount() == 0;
    }

    /**
     * Key function: tell this class the total number of jobs will provide data in the mapResultsQueue
     *
     * The total job count when we free threads blocked on waitForFinalReduce.  When we see numOfSubmittedJobs
     * MapResults from the queue, those threads are released.
     *
     * Until this function is called, those thread will block forever.  The numOfSubmittedJobs has a few constraints.
     * First, it must be >= 0.  0 indicates that in fact no jobs will ever be submitted (i.e., there's no
     * data coming) so the latch should be opened immediately.  If it's >= 1, we will wait until
     * we see numOfSubmittedJobs jobs before freeing them.
     *
     * Note that we throw an IllegalStateException if this function is called twice.
     *
     * @param numOfSubmittedJobs int >= 0 indicating the total number of MapResults that will
     *                           enqueue results into our queue
     */
    public synchronized void setTotalJobCount(final int numOfSubmittedJobs) {
        if ( numOfSubmittedJobs < 0 )
            throw new IllegalArgumentException("numOfSubmittedJobs must be >= 0, but saw " + numOfSubmittedJobs);
        if ( numJobsReduced > numOfSubmittedJobs )
            throw new IllegalArgumentException("numOfSubmittedJobs " + numOfSubmittedJobs + " < numJobsReduced " + numJobsReduced);
        if ( this.numSubmittedJobs != UNSET_NUM_SUBMITTED_JOBS)
            throw new IllegalStateException("setlastJobID called multiple times, but should only be called once");

        this.numSubmittedJobs = numOfSubmittedJobs;
        maybeReleaseLatch();
    }

    /**
     * Block until the last job has submitted its MapResult to our queue, and we've reduced it, and
     * return the reduce result resulting from applying reduce(...) to all MapResult elements.
     *
     * @return the total reduce result across all jobs
     * @throws InterruptedException
     */
    public ReduceType waitForFinalReduce() throws InterruptedException {
        countDownLatch.await();
        return sum;
    }
}
