package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.SimpleTimer;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * Reducer supporting two-threaded reduce of the map/reduce.
 *
 * The first thread, using the reduceAsMuchAsPossible function, actually reduces the data
 * as it arrives in the blockingQueue.
 *
 * The second thread, using the waitForFinalReduce, can block on this data structure
 * until that all jobs have arrived and been reduced.
 *
 * The key function for communication here is setLastJobID(), which the thread that submits
 * jobs that enqueue MapResults into the blocking queue must call ONCE to tell the
 * Reduce that ID of the last job that's been submitted.  When a job arrives with that
 * ID, this class frees a latch that allows thread blocked on waitForFinalReduce to proceed.
 *
 * This thread reads from mapResultsQueue until the poison EOF object arrives.  At each
 * stage is calls reduce(value, sum).  The blocking mapResultQueue ensures that the
 * queue waits until the mapResultQueue has a value to take. Then, it gets and waits
 * until the map result Future has a value.
 */
class Reducer<MapType, ReduceType> {
    private final static int UNSET_LAST_JOB_ID = -2;

    final CountDownLatch countDownLatch = new CountDownLatch(1);
    final NSReduceFunction<MapType, ReduceType> reduce;
    final SimpleTimer reduceTimer;

    /**
     * The sum of the reduce function applied to all MapResults.  After this Reducer
     * is done sum contains the final reduce result.
     */
    ReduceType sum;

    int lastJobID = UNSET_LAST_JOB_ID; // not yet set

    /**
     * The jobID of the last job we've seen
     */
    int prevJobID = -1; // no jobs observed

    /**
     * Create a new Reducer that will apply the reduce function with initialSum value
     * to values via reduceAsMuchAsPossible, timing the reduce function call costs with
     * reduceTimer
     *
     * @param reduce the reduce function to apply
     * @param reduceTimer the timer to time the reduce function call
     * @param initialSum the initial reduce sum
     */
    public Reducer(final NSReduceFunction<MapType, ReduceType> reduce,
                   final SimpleTimer reduceTimer,
                   final ReduceType initialSum) {
        if ( reduce == null ) throw new IllegalArgumentException("Reduce function cannot be null");
        if ( reduceTimer == null ) throw new IllegalArgumentException("reduceTimer cannot be null");

        this.reduce = reduce;
        this.reduceTimer = reduceTimer;
        this.sum = initialSum;
    }

    /**
     * Should we reduce the next value in the mapResultQueue?
     *
     *
     * @param mapResultQueue the queue of map results
     * @return true if we should reduce
     */
    @Requires("mapResultQueue != null")
    private synchronized boolean reduceNextValueInQueue(final PriorityBlockingQueue<MapResult<MapType>> mapResultQueue) {
        final MapResult<MapType> nextMapResult = mapResultQueue.peek();
        return nextMapResult != null && nextMapResult.getJobID() == prevJobID + 1;
    }

    /**
     * Reduce as much data as possible in mapResultQueue, returning the number of reduce calls completed
     *
     * As much as possible is defined as all of the MapResults in the queue are in order starting from the
     * lastJobID we reduced previously, up to the either the queue being empty or where the next MapResult
     * doesn't have JobID == prevJobID + 1.
     *
     * @param mapResultQueue a queue of MapResults in jobID order
     * @return the number of reduces run, from 0 >
     * @throws InterruptedException
     */
    @Ensures("result >= 0")
    public synchronized int reduceAsMuchAsPossible(final PriorityBlockingQueue<MapResult<MapType>> mapResultQueue) throws InterruptedException {
        if ( mapResultQueue == null ) throw new IllegalArgumentException("mapResultQueue cannot be null");
        int nReduces = 0;

        while ( reduceNextValueInQueue(mapResultQueue) ) {
            final MapResult<MapType> result = mapResultQueue.take();

            if ( result.getJobID() < prevJobID )
                // make sure the map results are coming in order
                throw new IllegalStateException("BUG: last jobID " + prevJobID + " > current jobID " + result.getJobID());

            prevJobID = result.getJobID();

            if ( ! result.isEOFMarker() ) {
                nReduces++;

                // apply reduce, keeping track of sum
                reduceTimer.restart();
                sum = reduce.apply(result.getValue(), sum);
                reduceTimer.stop();

            }

            maybeReleaseLatch();
        }

        return nReduces;
    }

    /**
     * release the latch if appropriate
     *
     * Appropriate means we've seen the last job, or there's only a single job id
     */
    private synchronized void maybeReleaseLatch() {
        if ( lastJobID != -2 && (prevJobID == lastJobID || lastJobID == -1) ) {
            // either we've already seen the last one prevJobID == lastJobID or
            // the last job ID is -1, meaning that no jobs were ever submitted
            countDownLatch.countDown();
        }
    }

    /**
     * For testing.
     * @return
     */
    protected synchronized boolean latchIsReleased() {
        return countDownLatch.getCount() == 0;
    }

    /**
     * Key function: tell this class the job ID of the last job that will provide data in the mapResultsQueue
     *
     * The last job id controls when we free threads blocked on waitForFinalReduce.  When we see the job
     * with this last job id, those threads are released.
     *
     * Until this function is called, those thread will block forever.  The last job id has a few constraints.
     * First, it must be >= -1.  -1 indicates that in fact no jobs will ever be submitted (i.e., there's no
     * data coming) so the latch should be opened immediately.  If it's >= 0, we will wait until
     * a job with that id arrives.
     *
     * Note that we throw an IllegalStateException if this function is called twice.
     *
     * @param lastJobID int >= -1 indicating the MapResult job id of the last job that will enqueue results into our queue
     */
    public synchronized void setLastJobID(final int lastJobID) {
        if ( lastJobID < -1 )
            throw new IllegalArgumentException("lastJobID must be > -1, but saw " + lastJobID);
        if ( this.lastJobID != UNSET_LAST_JOB_ID )
            throw new IllegalStateException("setlastJobID called multiple times, but should only be called once");

        this.lastJobID = lastJobID;
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
