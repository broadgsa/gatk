package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.AutoFormattingTime;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.threading.NamedThreadFactory;

import java.util.Iterator;
import java.util.List;
import java.util.concurrent.*;

/**
 * Framework for very fine grained MapReduce parallelism
 *
 * The overall framework works like this
 *
 * nano <- new Nanoschedule(inputBufferSize, numberOfMapElementsToProcessTogether, nThreads)
 * List[Input] outerData : outerDataLoop )
 *   result = nano.execute(outerData.iterator(), map, reduce)
 *
 * inputBufferSize determines how many elements from the input stream are read in one go by the
 * nanoscheduler.  The scheduler may hold up to inputBufferSize in memory at one time, as well
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
    private final static boolean LOG_MAP_TIMES = false;
    private final static boolean TIME_CALLS = true;

    private final static int MAP_BUFFER_SIZE_SCALE_FACTOR = 100;

    final int inputBufferSize;
    final int mapBufferSize;
    final int nThreads;
    final ExecutorService inputExecutor;
    final ExecutorService reduceExecutor;
    final ThreadPoolExecutor mapExecutor;

    boolean shutdown = false;
    boolean debug = false;
    private NSProgressFunction<InputType> progressFunction = null;

    final SimpleTimer outsideSchedulerTimer = TIME_CALLS ? new SimpleTimer("outside") : null;
    final SimpleTimer inputTimer = TIME_CALLS ? new SimpleTimer("input") : null;
    final SimpleTimer mapTimer = TIME_CALLS ? new SimpleTimer("map") : null;
    final SimpleTimer reduceTimer = TIME_CALLS ? new SimpleTimer("reduce") : null;

    /**
     * Create a new nanoscheduler with the desire characteristics requested by the argument
     *
     * @param inputBufferSize the number of input elements to read in each scheduling cycle.
     * @param nThreads the number of threads to use to get work done, in addition to the
     *                 thread calling execute
     */
    public NanoScheduler(final int inputBufferSize, final int nThreads) {
        if ( inputBufferSize < 1 ) throw new IllegalArgumentException("inputBufferSize must be >= 1, got " + inputBufferSize);
        if ( nThreads < 1 ) throw new IllegalArgumentException("nThreads must be >= 1, got " + nThreads);

        this.inputBufferSize = inputBufferSize;
        this.mapBufferSize = inputBufferSize * MAP_BUFFER_SIZE_SCALE_FACTOR;
        this.nThreads = nThreads;

        if ( nThreads == 1 ) {
            this.mapExecutor = null;
            this.inputExecutor = this.reduceExecutor = null;
        } else {
            this.mapExecutor = (ThreadPoolExecutor)Executors.newFixedThreadPool(nThreads-1, new NamedThreadFactory("NS-map-thread-%d"));
            this.mapExecutor.setRejectedExecutionHandler(new ThreadPoolExecutor.CallerRunsPolicy());
            this.inputExecutor = Executors.newSingleThreadExecutor(new NamedThreadFactory("NS-input-thread-%d"));
            this.reduceExecutor = Executors.newSingleThreadExecutor(new NamedThreadFactory("NS-reduce-thread-%d"));
        }

        // start timing the time spent outside of the nanoScheduler
        outsideSchedulerTimer.start();
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
    public int getInputBufferSize() {
        return inputBufferSize;
    }

    /**
     * Tells this nanoScheduler to shutdown immediately, releasing all its resources.
     *
     * After this call, execute cannot be invoked without throwing an error
     */
    public void shutdown() {
        outsideSchedulerTimer.stop();

        if ( nThreads > 1 ) {
            shutdownExecutor("inputExecutor", inputExecutor);
            shutdownExecutor("mapExecutor", mapExecutor);
            shutdownExecutor("reduceExecutor", reduceExecutor);
        }
        shutdown = true;

        if (TIME_CALLS) {
            printTimerInfo("Input   time", inputTimer);
            printTimerInfo("Map     time", mapTimer);
            printTimerInfo("Reduce  time", reduceTimer);
            printTimerInfo("Outside time", outsideSchedulerTimer);
        }
    }

    /**
     * Helper function to cleanly shutdown an execution service, checking that the execution
     * state is clean when it's done.
     *
     * @param name a string name for error messages for the executorService we are shutting down
     * @param executorService the executorService to shut down
     */
    @Requires({"name != null", "executorService != null"})
    @Ensures("executorService.isShutdown()")
    private void shutdownExecutor(final String name, final ExecutorService executorService) {
        if ( executorService.isShutdown() || executorService.isTerminated() )
            throw new IllegalStateException("Executor service " + name + " is already shut down!");

        final List<Runnable> remaining = executorService.shutdownNow();
        if ( ! remaining.isEmpty() )
            throw new IllegalStateException(remaining.size() + " remaining tasks found in an executor " + name + ", unexpected behavior!");
    }

    /**
     * Print to logger.info timing information from timer, with name label
     *
     * @param label the name of the timer to display.  Should be human readable
     * @param timer the timer whose elapsed time we will display
     */
    @Requires({"label != null", "timer != null"})
    private void printTimerInfo(final String label, final SimpleTimer timer) {
        final double total = inputTimer.getElapsedTime() + mapTimer.getElapsedTime()
                + reduceTimer.getElapsedTime() + outsideSchedulerTimer.getElapsedTime();
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

    /**
     * @return are we displaying verbose debugging information about the scheduling?
     */
    public boolean isDebug() {
        return debug;
    }

    /**
     * Helper function to display a String.formatted message if we are doing verbose debugging
     *
     * @param format the format argument suitable for String.format
     * @param args the arguments for String.format
     */
    @Requires("format != null")
    private void debugPrint(final String format, Object ... args) {
        if ( isDebug() )
            logger.info("Thread " + Thread.currentThread().getId() + ":" + String.format(format, args));
    }

    /**
     * Turn on/off verbose debugging
     *
     * @param debug true if we want verbose debugging
     */
    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    /**
     * Set the progress callback function to progressFunction
     *
     * The progress callback is invoked after each buffer size elements have been processed by map/reduce
     *
     * @param progressFunction a progress function to call, or null if you don't want any progress callback
     */
    public void setProgressFunction(final NSProgressFunction<InputType> progressFunction) {
        this.progressFunction = progressFunction;
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
     * Note that this function goes through a single threaded fast path if the number of threads
     * is 1.
     *
     * @param inputReader an iterator providing us with the input data to nanoSchedule map/reduce over
     * @param map the map function from input type -> map type, will be applied in parallel to each input
     * @param reduce the reduce function from map type + reduce type -> reduce type to be applied in order to map results
     * @return the last reduce value
     */
    public ReduceType execute(final Iterator<InputType> inputReader,
                              final NSMapFunction<InputType, MapType> map,
                              final ReduceType initialValue,
                              final NSReduceFunction<MapType, ReduceType> reduce) {
        if ( isShutdown() ) throw new IllegalStateException("execute called on already shutdown NanoScheduler");
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( map == null ) throw new IllegalArgumentException("map function cannot be null");
        if ( reduce == null ) throw new IllegalArgumentException("reduce function cannot be null");

        outsideSchedulerTimer.stop();

        ReduceType result;
        if ( ALLOW_SINGLE_THREAD_FASTPATH && getnThreads() == 1 ) {
            result = executeSingleThreaded(inputReader, map, initialValue, reduce);
        } else {
            result = executeMultiThreaded(inputReader, map, initialValue, reduce);
        }

        outsideSchedulerTimer.restart();
        return result;
    }

    /**
     * Simple efficient reference implementation for single threaded execution.
     *
     * @return the reduce result of this map/reduce job
     */
    @Requires({"inputReader != null", "map != null", "reduce != null"})
    private ReduceType executeSingleThreaded(final Iterator<InputType> inputReader,
                                             final NSMapFunction<InputType, MapType> map,
                                             final ReduceType initialValue,
                                             final NSReduceFunction<MapType, ReduceType> reduce) {
        ReduceType sum = initialValue;
        int i = 0;

        // start timer to ensure that both hasNext and next are caught by the timer
        if ( TIME_CALLS ) inputTimer.restart();
        while ( inputReader.hasNext() ) {
            final InputType input = inputReader.next();
            if ( TIME_CALLS ) inputTimer.stop();

            // map
            if ( TIME_CALLS ) mapTimer.restart();
            final long preMapTime = LOG_MAP_TIMES ? 0 : mapTimer.currentTimeNano();
            final MapType mapValue = map.apply(input);
            if ( LOG_MAP_TIMES ) logger.info("MAP TIME " + (mapTimer.currentTimeNano() - preMapTime));
            if ( TIME_CALLS ) mapTimer.stop();

            if ( i++ % inputBufferSize == 0 && progressFunction != null )
                progressFunction.progress(input);

            // reduce
            if ( TIME_CALLS ) reduceTimer.restart();
            sum = reduce.apply(mapValue, sum);
            if ( TIME_CALLS ) reduceTimer.stop();

            if ( TIME_CALLS ) inputTimer.restart();
        }

        return sum;
    }

    /**
     * Efficient parallel version of Map/Reduce
     *
     * @return the reduce result of this map/reduce job
     */
    @Requires({"inputReader != null", "map != null", "reduce != null"})
    private ReduceType executeMultiThreaded(final Iterator<InputType> inputReader,
                                            final NSMapFunction<InputType, MapType> map,
                                            final ReduceType initialValue,
                                            final NSReduceFunction<MapType, ReduceType> reduce) {
        debugPrint("Executing nanoScheduler");

        // a blocking queue that limits the number of input datum to the requested buffer size
        final BlockingQueue<InputProducer<InputType>.InputValue> inputQueue
                = new LinkedBlockingDeque<InputProducer<InputType>.InputValue>(inputBufferSize);

        // a priority queue that stores up to mapBufferSize elements
        // produced by completed map jobs.
        final BlockingQueue<Future<MapResult<MapType>>> mapResultQueue =
                new LinkedBlockingDeque<Future<MapResult<MapType>>>(mapBufferSize);

        // Start running the input reader thread
        inputExecutor.submit(new InputProducer<InputType>(inputReader, inputTimer, inputQueue));

        // Start running the reducer thread
        final ReducerThread<MapType, ReduceType> reducer
                = new ReducerThread<MapType, ReduceType>(reduce, reduceTimer, initialValue, mapResultQueue);
        final Future<ReduceType> reduceResult = reduceExecutor.submit(reducer);

        try {
            int numJobs = 0;

            while ( true ) {
                // block on input
                final InputProducer<InputType>.InputValue inputEnqueueWrapped = inputQueue.take();

                if ( ! inputEnqueueWrapped.isLast() ) {
                    // get the object itself
                    final InputType input = inputEnqueueWrapped.getValue();

                    // the next map call has jobID + 1
                    numJobs++;

                    // send job for map via the completion service
                    final CallableMap doMap = new CallableMap(map, numJobs, input);
                    final Future<MapResult<MapType>> mapJob = mapExecutor.submit(doMap);
                    mapResultQueue.put(mapJob);

                    debugPrint("  Done with cycle of map/reduce");

                    if ( numJobs % inputBufferSize == 0 && progressFunction != null )
                        progressFunction.progress(input);
                } else {
                    mapResultQueue.put(new FutureValue<MapResult<MapType>>(new MapResult<MapType>()));
                    return reduceResult.get(); // wait for our result of reduce
                }
            }
        } catch (InterruptedException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        } catch (ExecutionException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        }
    }

    /**
     * A simple callable version of the map function for use with the executor pool
     */
    private class CallableMap implements Callable<MapResult<MapType>> {
        final int id;
        final InputType input;
        final NSMapFunction<InputType, MapType> map;

        @Requires({"map != null"})
        private CallableMap(final NSMapFunction<InputType, MapType> map,
                            final int id,
                            final InputType input) {
            this.id = id;
            this.input = input;
            this.map = map;
        }

        @Override
        public MapResult<MapType> call() {
            if ( TIME_CALLS ) mapTimer.restart();
            if ( debug ) debugPrint("\t\tmap " + input);
            final MapType result = map.apply(input);
            if ( TIME_CALLS ) mapTimer.stop();
            return new MapResult<MapType>(result, id);
        }
    }
}
