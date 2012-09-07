package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Thread that runs the reduce of the map/reduce.
 *
 * This thread reads from mapResultsQueue until the poison EOF object arrives.  At each
 * stage is calls reduce(value, sum).  The blocking mapResultQueue ensures that the
 * queue waits until the mapResultQueue has a value to take. Then, it gets and waits
 * until the map result Future has a value.
 */
class ReducerThread<MapType, ReduceType> implements Callable<ReduceType> {
    final NSReduceFunction<MapType, ReduceType> reduce;
    final SimpleTimer reduceTimer;
    final BlockingQueue<Future<MapResult<MapType>>> mapResultQueue;

    ReduceType sum;
    int lastJobID = -1;

    public ReducerThread(final NSReduceFunction<MapType, ReduceType> reduce,
                         final SimpleTimer reduceTimer,
                         final ReduceType sum,
                         final BlockingQueue<Future<MapResult<MapType>>> mapResultQueue) {
        if ( reduce == null ) throw new IllegalArgumentException("Reduce function cannot be null");
        if ( mapResultQueue == null ) throw new IllegalArgumentException("mapResultQueue cannot be null");

        this.reduce = reduce;
        this.reduceTimer = reduceTimer;
        this.sum = sum;
        this.mapResultQueue = mapResultQueue;
    }

    public ReduceType call() {
        try {
            while ( true ) {
                final MapResult<MapType> result = mapResultQueue.take().get();
                if ( result.isLast() ) {
                    // we are done, just return sum
                    return sum;
                }
                else if ( result.getJobID() < lastJobID ) {
                    // make sure the map results are coming in order
                    throw new IllegalStateException("BUG: last jobID " + lastJobID + " > current jobID " + result.getJobID());
                } else {
                    // apply reduce, keeping track of sum
                    if ( reduceTimer != null ) reduceTimer.restart();
                    sum = reduce.apply(result.getValue(), sum);
                    if ( reduceTimer != null ) reduceTimer.stop();
                }
            }
        } catch (ExecutionException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        } catch (InterruptedException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        }
    }
}
