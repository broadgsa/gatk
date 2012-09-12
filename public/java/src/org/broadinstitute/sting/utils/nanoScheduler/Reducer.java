package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.utils.SimpleTimer;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Thread that runs the reduce of the map/reduce.
 *
 * This thread reads from mapResultsQueue until the poison EOF object arrives.  At each
 * stage is calls reduce(value, sum).  The blocking mapResultQueue ensures that the
 * queue waits until the mapResultQueue has a value to take. Then, it gets and waits
 * until the map result Future has a value.
 */
class Reducer<MapType, ReduceType> {
    final CountDownLatch countDownLatch = new CountDownLatch(1);

    final NSReduceFunction<MapType, ReduceType> reduce;
    final SimpleTimer reduceTimer;

    ReduceType sum;
    int lastJobID = -2; // not yet set
    int prevJobID = -1; // no jobs observed

    public Reducer(final NSReduceFunction<MapType, ReduceType> reduce,
                   final SimpleTimer reduceTimer,
                   final ReduceType initialSum) {
        if ( reduce == null ) throw new IllegalArgumentException("Reduce function cannot be null");
        if ( reduceTimer == null ) throw new IllegalArgumentException("reduceTimer cannot be null");

        this.reduce = reduce;
        this.reduceTimer = reduceTimer;
        this.sum = initialSum;
    }

    private synchronized boolean readyToReduce(final BlockingQueue<MapResult<MapType>> mapResultQueue) {
        final MapResult<MapType> nextMapResult = mapResultQueue.peek();
        return nextMapResult != null && nextMapResult.getJobID() == prevJobID + 1;
    }

    public synchronized int reduceAsMuchAsPossible(final BlockingQueue<MapResult<MapType>> mapResultQueue) throws InterruptedException {
        int nReduces = 0;

        while ( readyToReduce(mapResultQueue) ) {
            final MapResult<MapType> result = mapResultQueue.take();

            if ( result.getJobID() < prevJobID )
                // make sure the map results are coming in order
                throw new IllegalStateException("BUG: last jobID " + prevJobID + " > current jobID " + result.getJobID());

            prevJobID = result.getJobID();

            if ( ! result.isLast() ) { // TODO -- rename to isEmpty
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

    private synchronized void maybeReleaseLatch() {
        if ( lastJobID != -2 && (prevJobID == lastJobID || lastJobID == -1) ) {
            // either we've already seen the last one prevJobID == lastJobID or
            // the last job ID is -1, meaning that no jobs were ever submitted
            countDownLatch.countDown();
        }
    }

    public synchronized void setLastJobID(final int lastJobID) {
        if ( lastJobID < -1 ) throw new IllegalArgumentException("lastJobID must be > -1, but saw " + lastJobID);
        this.lastJobID = lastJobID;
        maybeReleaseLatch();
    }

    public ReduceType waitForFinalReduce() throws InterruptedException {
        countDownLatch.await();
        return sum;
    }
}
