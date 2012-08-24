package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.*;

/**
 * Framework for very fine grained MapReduce parallelism
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:47 AM
 */
public class NanoScheduler<InputType, MapType, ReduceType> {
    final int bufferSize;
    final int nThreads;
    final Iterator<InputType> inputReader;
    final MapFunction<InputType, MapType> map;
    final ReduceFunction<MapType, ReduceType> reduce;

    public NanoScheduler(final int bufferSize,
                         final int nThreads,
                         final Iterator<InputType> inputReader,
                         final MapFunction<InputType, MapType> map,
                         final ReduceFunction<MapType, ReduceType> reduce) {
        if ( bufferSize < 1 ) throw new IllegalArgumentException("bufferSize must be >= 1, got " + bufferSize);
        if ( nThreads < 1 ) throw new IllegalArgumentException("nThreads must be >= 1, got " + nThreads);

        this.bufferSize = bufferSize;
        this.inputReader = inputReader;
        this.map = map;
        this.reduce = reduce;
        this.nThreads = nThreads;
    }

    public int getnThreads() {
        return nThreads;
    }

    private int getBufferSize() {
        return bufferSize;
    }

    public ReduceType execute() {
        if ( getnThreads() == 1 ) {
            return executeSingleThreaded();
        } else {
            return executeMultiThreaded();
        }
    }

    /**
     * Simple efficient reference implementation for single threaded execution
     * @return the reduce result of this map/reduce job
     */
    private ReduceType executeSingleThreaded() {
        ReduceType sum = reduce.init();
        while ( inputReader.hasNext() ) {
            final InputType input = inputReader.next();
            final MapType mapValue = map.apply(input);
            sum = reduce.apply(mapValue, sum);
        }
        return sum;
    }

    /**
     * Efficient parallel version of Map/Reduce
     *
     * @return the reduce result of this map/reduce job
     */
    private ReduceType executeMultiThreaded() {
        final ExecutorService executor = Executors.newFixedThreadPool(getnThreads() - 1);

        ReduceType sum = reduce.init();
        while ( inputReader.hasNext() ) {
            try {
                // read in our input values
                final Queue<InputType> inputs = readInputs();

                // send jobs for map
                final Queue<Future<MapType>> mapQueue = submitMapJobs(executor, inputs);

                // send off the reduce job, and block until we get at least one reduce result
                sum = reduceParallel(mapQueue, sum);
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            } catch (ExecutionException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            }
        }

        final List<Runnable> remaining = executor.shutdownNow();
        if ( ! remaining.isEmpty() )
            throw new ReviewedStingException("Remaining tasks found in the executor, unexpected behavior!");

        return sum;
    }

    @Requires("! mapQueue.isEmpty()")
    private ReduceType reduceParallel(final Queue<Future<MapType>> mapQueue, final ReduceType initSum)
            throws InterruptedException, ExecutionException {
        ReduceType sum = initSum;

        // while mapQueue has something in it to reduce
        for ( final Future<MapType> future : mapQueue ) {
            // block until we get the value for this task
            final MapType value = future.get();
            sum = reduce.apply(value, sum);
        }

        return sum;
    }

    /**
     * Read up to inputBufferSize elements from inputReader
     *
     * @return a queue of inputs read in, containing one or more values of InputType read in
     */
    @Requires("inputReader.hasNext()")
    @Ensures("!result.isEmpty()")
    private Queue<InputType> readInputs() {
        int n = 0;
        final Queue<InputType> inputs = new LinkedList<InputType>();
        while ( inputReader.hasNext() && n < getBufferSize() ) {
            final InputType input = inputReader.next();
            inputs.add(input);
            n++;
        }
        return inputs;
    }

    @Ensures("result.size() == inputs.size()")
    private Queue<Future<MapType>> submitMapJobs(final ExecutorService executor, final Queue<InputType> inputs) {
        final Queue<Future<MapType>> mapQueue = new LinkedList<Future<MapType>>();

        for ( final InputType input : inputs ) {
            final CallableMap doMap = new CallableMap(input);
            final Future<MapType> future = executor.submit(doMap);
            mapQueue.add(future);
        }

        return mapQueue;
    }

    /**
     * A simple callable version of the map function for use with the executor pool
     */
    private class CallableMap implements Callable<MapType> {
        final InputType input;

        private CallableMap(final InputType input) {
            this.input = input;
        }

        @Override public MapType call() throws Exception {
            return map.apply(input);
        }
    }
}
