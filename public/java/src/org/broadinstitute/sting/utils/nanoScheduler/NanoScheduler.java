package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.*;

/**
 * Framework for very fine grained MapReduce parallelism
 *
 * The overall framework works like this
 *
 * nano <- new Nanoschedule(bufferSize, numberOfMapElementsToProcessTogether, nThreads)
 * List[Input] outerData : outerDataLoop )
 *   result = nano.execute(outerData.iterator(), map, reduce)
 *
 * bufferSize determines how many elements from the input stream are read in one go by the
 * nanoscheduler.  The scheduler may hold up to bufferSize in memory at one time, as well
 * as up to inputBufferSize map results as well.
 *
 * numberOfMapElementsToProcessTogether determines how many input elements are processed
 * together each thread cycle.  For example, if this value is 10, then the input data
 * is grouped together in units of 10 elements each, and map called on each in term.  The more
 * heavy-weight the map function is, in terms of CPU costs, the more it makes sense to
 * have this number be small.  The lighter the CPU cost per element, though, the more this
 * parameter introduces overhead due to need to context switch among threads to process
 * each input element.  A value of -1 lets the nanoscheduler guess at a reasonable trade-off value.
 *
 * nThreads is a bit obvious yes?  Note though that the nanoscheduler assumes that it gets 1 thread
 * from its client during the execute call, as this call blocks until all work is done.  The caller
 * thread is put to work by execute to help with the processing of the data.  So in reality the
 * nanoScheduler only spawn nThreads - 1 additional workers (if this is > 1).
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:47 AM
 */
public class NanoScheduler<InputType, MapType, ReduceType> {
    private static Logger logger = Logger.getLogger(NanoScheduler.class);

    final int bufferSize;
    final int mapGroupSize;
    final int nThreads;
    final ExecutorService executor;
    boolean shutdown = false;
    boolean debug = false;

    /**
     * Create a new nanoschedule with the desire characteristics requested by the argument
     *
     * @param bufferSize the number of input elements to read in each scheduling cycle.
     * @param mapGroupSize How many inputs should be grouped together per map?  If -1 we make a reasonable guess
     * @param nThreads the number of threads to use to get work done, in addition to the thread calling execute
     */
    public NanoScheduler(final int bufferSize,
                         final int mapGroupSize,
                         final int nThreads) {
        if ( bufferSize < 1 ) throw new IllegalArgumentException("bufferSize must be >= 1, got " + bufferSize);
        if ( nThreads < 1 ) throw new IllegalArgumentException("nThreads must be >= 1, got " + nThreads);

        if ( mapGroupSize > bufferSize ) throw new IllegalArgumentException("mapGroupSize " + mapGroupSize + " must be <= bufferSize " + bufferSize);
        if ( mapGroupSize == 0 || mapGroupSize < -1 ) throw new IllegalArgumentException("mapGroupSize cannot be <= 0" + mapGroupSize);

        this.bufferSize = bufferSize;
        this.nThreads = nThreads;

        if ( mapGroupSize == -1 ) {
            this.mapGroupSize = (int)Math.ceil(this.bufferSize / (10.0*this.nThreads));
            logger.info(String.format("Dynamically setting grouping size to %d based on buffer size %d and n threads %d",
                    this.mapGroupSize, this.bufferSize, this.nThreads));
        } else {
            this.mapGroupSize = mapGroupSize;
        }

        this.executor = nThreads == 1 ? null : Executors.newFixedThreadPool(nThreads - 1);
    }

    /**
     * The number of parallel map threads in use with this NanoScheduler
     * @return
     */
    @Ensures("result > 0")
    public int getnThreads() {
        return nThreads;
    }

    /**
     * The input buffer size used by this NanoScheduler
     * @return
     */
    @Ensures("result > 0")
    public int getBufferSize() {
        return bufferSize;
    }

    /**
     * The grouping size used by this NanoScheduler
     * @return
     */
    @Ensures("result > 0")
    public int getMapGroupSize() {
        return mapGroupSize;
    }

    /**
     * Tells this nanoScheduler to shutdown immediately, releasing all its resources.
     *
     * After this call, execute cannot be invoked without throwing an error
     */
    public void shutdown() {
        if ( executor != null ) {
            final List<Runnable> remaining = executor.shutdownNow();
            if ( ! remaining.isEmpty() )
                throw new IllegalStateException("Remaining tasks found in the executor, unexpected behavior!");
        }
        shutdown = true;
    }

    /**
     * @return true if this nanoScheduler is shutdown, or false if its still open for business
     */
    public boolean isShutdown() {
        return shutdown;
    }

    public boolean isDebug() {
        return debug;
    }

    private void debugPrint(final String format, Object ... args) {
        if ( isDebug() )
            logger.info("Thread " + Thread.currentThread().getId() + ":" + String.format(format, args));
    }


    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    /**
     * Execute a map/reduce job with this nanoScheduler
     *
     * Data comes from inputReader.  Will be read until hasNext() == false.
     * map is called on each element provided by inputReader.  No order of operations is guarenteed
     * reduce is called in order of the input data provided by inputReader on the result of map() applied
     * to each element.
     *
     * Note that the caller thread is put to work with this function call.  The call doesn't return
     * until all elements have been processes.
     *
     * It is safe to call this function repeatedly on a single nanoScheduler, at least until the
     * shutdown method is called.
     *
     * @param inputReader
     * @param map
     * @param reduce
     * @return
     */
    public ReduceType execute(final Iterator<InputType> inputReader,
                              final MapFunction<InputType, MapType> map,
                              final ReduceType initialValue,
                              final ReduceFunction<MapType, ReduceType> reduce) {
        if ( isShutdown() ) throw new IllegalStateException("execute called on already shutdown NanoScheduler");
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( map == null ) throw new IllegalArgumentException("map function cannot be null");
        if ( reduce == null ) throw new IllegalArgumentException("reduce function cannot be null");

        if ( getnThreads() == 1 ) {
            return executeSingleThreaded(inputReader, map, initialValue, reduce);
        } else {
            return executeMultiThreaded(inputReader, map, initialValue, reduce);
        }
    }

    /**
     * Simple efficient reference implementation for single threaded execution
     * @return the reduce result of this map/reduce job
     */
    private ReduceType executeSingleThreaded(final Iterator<InputType> inputReader,
                                             final MapFunction<InputType, MapType> map,
                                             final ReduceType initialValue,
                                             final ReduceFunction<MapType, ReduceType> reduce) {
        ReduceType sum = initialValue;
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
    private ReduceType executeMultiThreaded(final Iterator<InputType> inputReader,
                                            final MapFunction<InputType, MapType> map,
                                            final ReduceType initialValue,
                                            final ReduceFunction<MapType, ReduceType> reduce) {
        debugPrint("Executing nanoScheduler");
        ReduceType sum = initialValue;
        while ( inputReader.hasNext() ) {
            try {
                // read in our input values
                final List<InputType> inputs = readInputs(inputReader);

                // send jobs for map
                final Queue<Future<List<MapType>>> mapQueue = submitMapJobs(map, executor, inputs);

                // send off the reduce job, and block until we get at least one reduce result
                sum = reduceParallel(reduce, mapQueue, sum);
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            } catch (ExecutionException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            }
        }

        return sum;
    }

    @Requires({"reduce != null", "! mapQueue.isEmpty()"})
    private ReduceType reduceParallel(final ReduceFunction<MapType, ReduceType> reduce,
                                      final Queue<Future<List<MapType>>> mapQueue,
                                      final ReduceType initSum)
            throws InterruptedException, ExecutionException {
        ReduceType sum = initSum;

        // while mapQueue has something in it to reduce
        for ( final Future<List<MapType>> future : mapQueue ) {
            for ( final MapType value : future.get() ) // block until we get the values for this task
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
    private List<InputType> readInputs(final Iterator<InputType> inputReader) {
        int n = 0;
        final List<InputType> inputs = new LinkedList<InputType>();
        while ( inputReader.hasNext() && n < getBufferSize() ) {
            final InputType input = inputReader.next();
            inputs.add(input);
            n++;
        }
        return inputs;
    }

    @Requires({"map != null", "! inputs.isEmpty()"})
    private Queue<Future<List<MapType>>> submitMapJobs(final MapFunction<InputType, MapType> map,
                                                       final ExecutorService executor,
                                                       final List<InputType> inputs) {
        final Queue<Future<List<MapType>>> mapQueue = new LinkedList<Future<List<MapType>>>();

        for ( final List<InputType> subinputs : Utils.groupList(inputs, getMapGroupSize()) ) {
            final CallableMap doMap = new CallableMap(map, subinputs);
            final Future<List<MapType>> future = executor.submit(doMap);
            mapQueue.add(future);
        }

        return mapQueue;
    }

    /**
     * A simple callable version of the map function for use with the executor pool
     */
    private class CallableMap implements Callable<List<MapType>> {
        final List<InputType> inputs;
        final MapFunction<InputType, MapType> map;

        @Requires({"map != null", "inputs.size() <= getMapGroupSize()"})
        private CallableMap(final MapFunction<InputType, MapType> map, final List<InputType> inputs) {
            this.inputs = inputs;
            this.map = map;
        }

        @Ensures("result.size() == inputs.size()")
        @Override public List<MapType> call() throws Exception {
            final List<MapType> outputs = new LinkedList<MapType>();
            for ( final InputType input : inputs )
                outputs.add(map.apply(input));
            debugPrint("    Processed %d elements with map", outputs.size());
            return outputs;
        }
    }
}
