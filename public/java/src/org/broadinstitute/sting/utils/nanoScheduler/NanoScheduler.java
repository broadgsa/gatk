package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.AutoFormattingTime;
import org.broadinstitute.sting.utils.SimpleTimer;
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
    private final static Logger logger = Logger.getLogger(NanoScheduler.class);
    private final static boolean ALLOW_SINGLE_THREAD_FASTPATH = true;
    private final static boolean TIME_CALLS = true;

    final int bufferSize;
    final int nThreads;
    final ExecutorService executor;
    boolean shutdown = false;
    boolean debug = false;

    final SimpleTimer inputTimer = new SimpleTimer();
    final SimpleTimer mapTimer = new SimpleTimer();
    final SimpleTimer reduceTimer = new SimpleTimer();

    /**
     * Create a new nanoschedule with the desire characteristics requested by the argument
     *
     * @param bufferSize the number of input elements to read in each scheduling cycle.
     * @param nThreads the number of threads to use to get work done, in addition to the thread calling execute
     */
    public NanoScheduler(final int bufferSize,
                         final int nThreads) {
        if ( bufferSize < 1 ) throw new IllegalArgumentException("bufferSize must be >= 1, got " + bufferSize);
        if ( nThreads < 1 ) throw new IllegalArgumentException("nThreads must be >= 1, got " + nThreads);

        this.bufferSize = bufferSize;
        this.nThreads = nThreads;
        this.executor = nThreads == 1 ? null : Executors.newFixedThreadPool(nThreads);
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

        if (TIME_CALLS) {
            printTimerInfo("Input  time", inputTimer);
            printTimerInfo("Map    time", mapTimer);
            printTimerInfo("Reduce time", reduceTimer);
        }
    }

    private void printTimerInfo(final String label, final SimpleTimer timer) {
        final double total = inputTimer.getElapsedTime() + mapTimer.getElapsedTime() + reduceTimer.getElapsedTime();
        final double myTimeInSec = timer.getElapsedTime();
        final double myTimePercent = myTimeInSec / total * 100;
        logger.info(String.format("%s: %s (%5.2f%%)", label, new AutoFormattingTime(myTimeInSec), myTimePercent));
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
     * @param inputReader an iterator providing us with the input data to nanoSchedule map/reduce over
     * @param map the map function from input type -> map type, will be applied in parallel to each input
     * @param reduce the reduce function from map type + reduce type -> reduce type to be applied in order to map results
     * @return the last reduce value
     */
    public ReduceType execute(final Iterator<InputType> inputReader,
                              final MapFunction<InputType, MapType> map,
                              final ReduceType initialValue,
                              final ReduceFunction<MapType, ReduceType> reduce) {
        if ( isShutdown() ) throw new IllegalStateException("execute called on already shutdown NanoScheduler");
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( map == null ) throw new IllegalArgumentException("map function cannot be null");
        if ( reduce == null ) throw new IllegalArgumentException("reduce function cannot be null");

        if ( ALLOW_SINGLE_THREAD_FASTPATH && getnThreads() == 1 ) {
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
                final Queue<Future<MapType>> mapQueue = submitMapJobs(map, executor, inputs);

                // send off the reduce job, and block until we get at least one reduce result
                sum = reduceSerial(reduce, mapQueue, sum);
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            } catch (ExecutionException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            }
        }

        return sum;
    }

    @Requires({"reduce != null", "! mapQueue.isEmpty()"})
    private ReduceType reduceSerial(final ReduceFunction<MapType, ReduceType> reduce,
                                    final Queue<Future<MapType>> mapQueue,
                                    final ReduceType initSum)
            throws InterruptedException, ExecutionException {
        ReduceType sum = initSum;

        // while mapQueue has something in it to reduce
        for ( final Future<MapType> future : mapQueue ) {
            final MapType value = future.get(); // block until we get the values for this task

            if ( TIME_CALLS) reduceTimer.restart();
            sum = reduce.apply(value, sum);
            if ( TIME_CALLS) reduceTimer.stop();
        }

        return sum;
    }

    /**
     * Read up to inputBufferSize elements from inputReader
     *
     * @return a queue of input read in, containing one or more values of InputType read in
     */
    @Requires("inputReader.hasNext()")
    @Ensures("!result.isEmpty()")
    private List<InputType> readInputs(final Iterator<InputType> inputReader) {
        int n = 0;
        final List<InputType> inputs = new LinkedList<InputType>();

        if ( TIME_CALLS) inputTimer.restart();
        while ( inputReader.hasNext() && n < getBufferSize() ) {
            final InputType input = inputReader.next();
            inputs.add(input);
            n++;
        }
        if ( TIME_CALLS) inputTimer.stop();

        return inputs;
    }

    @Requires({"map != null", "! inputs.isEmpty()"})
    private Queue<Future<MapType>> submitMapJobs(final MapFunction<InputType, MapType> map,
                                                 final ExecutorService executor,
                                                 final List<InputType> inputs) {
        final Queue<Future<MapType>> mapQueue = new LinkedList<Future<MapType>>();

        for ( final InputType input : inputs ) {
            final CallableMap doMap = new CallableMap(map, input);
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
        final MapFunction<InputType, MapType> map;

        @Requires({"map != null"})
        private CallableMap(final MapFunction<InputType, MapType> map, final InputType inputs) {
            this.input = inputs;
            this.map = map;
        }

        @Override public MapType call() throws Exception {
            if ( TIME_CALLS) mapTimer.restart();
            final MapType result = map.apply(input);
            if ( TIME_CALLS) mapTimer.stop();
            return result;
        }
    }
}
