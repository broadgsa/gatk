package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.AutoFormattingTime;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.collections.Pair;
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
    private final static boolean LOG_MAP_TIMES = false;
    private final static boolean TIME_CALLS = true;

    final int bufferSize;
    final int nThreads;
    final ExecutorService inputExecutor;
    final ExecutorService mapExecutor;
    boolean shutdown = false;
    boolean debug = false;

    private NanoSchedulerProgressFunction<InputType> progressFunction = null;

    final SimpleTimer outsideSchedulerTimer = new SimpleTimer("outside");
    final SimpleTimer inputTimer = new SimpleTimer("input");
    final SimpleTimer mapTimer = new SimpleTimer("map");
    final SimpleTimer reduceTimer = new SimpleTimer("reduce");

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
        this.mapExecutor = nThreads == 1 ? null : Executors.newFixedThreadPool(nThreads-1);
        this.inputExecutor = Executors.newSingleThreadExecutor();

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
    public int getBufferSize() {
        return bufferSize;
    }

    /**
     * Tells this nanoScheduler to shutdown immediately, releasing all its resources.
     *
     * After this call, execute cannot be invoked without throwing an error
     */
    public void shutdown() {
        outsideSchedulerTimer.stop();

        if ( mapExecutor != null ) {
            final List<Runnable> remaining = mapExecutor.shutdownNow();
            if ( ! remaining.isEmpty() )
                throw new IllegalStateException("Remaining tasks found in the mapExecutor, unexpected behavior!");
        }
        shutdown = true;

        if (TIME_CALLS) {
            printTimerInfo("Input   time", inputTimer);
            printTimerInfo("Map     time", mapTimer);
            printTimerInfo("Reduce  time", reduceTimer);
            printTimerInfo("Outside time", outsideSchedulerTimer);
        }
    }

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
     * Set the progress callback function to progressFunction
     *
     * The progress callback is invoked after each buffer size elements have been processed by map/reduce
     *
     * @param progressFunction a progress function to call, or null if you don't want any progress callback
     */
    public void setProgressFunction(final NanoSchedulerProgressFunction<InputType> progressFunction) {
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
     * @param inputReader an iterator providing us with the input data to nanoSchedule map/reduce over
     * @param map the map function from input type -> map type, will be applied in parallel to each input
     * @param reduce the reduce function from map type + reduce type -> reduce type to be applied in order to map results
     * @return the last reduce value
     */
    public ReduceType execute(final Iterator<InputType> inputReader,
                              final NanoSchedulerMapFunction<InputType, MapType> map,
                              final ReduceType initialValue,
                              final NanoSchedulerReduceFunction<MapType, ReduceType> reduce) {
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
     * Simple efficient reference implementation for single threaded execution
     * @return the reduce result of this map/reduce job
     */
    private ReduceType executeSingleThreaded(final Iterator<InputType> inputReader,
                                             final NanoSchedulerMapFunction<InputType, MapType> map,
                                             final ReduceType initialValue,
                                             final NanoSchedulerReduceFunction<MapType, ReduceType> reduce) {
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

            if ( i++ % bufferSize == 0 && progressFunction != null )
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
    private ReduceType executeMultiThreaded(final Iterator<InputType> inputReader,
                                            final NanoSchedulerMapFunction<InputType, MapType> map,
                                            final ReduceType initialValue,
                                            final NanoSchedulerReduceFunction<MapType, ReduceType> reduce) {
        debugPrint("Executing nanoScheduler");
        ReduceType sum = initialValue;
        boolean done = false;

        final BlockingQueue<InputDatum> inputQueue = new LinkedBlockingDeque<InputDatum>(bufferSize);

        inputExecutor.submit(new InputProducer(inputReader, inputQueue));

        while ( ! done ) {
            try {
                final Pair<List<InputType>, Boolean> readResults = readInputs(inputQueue);
                final List<InputType> inputs = readResults.getFirst();
                done = readResults.getSecond();

                if ( ! inputs.isEmpty() ) {
                    // send jobs for map
                    final Queue<Future<MapType>> mapQueue = submitMapJobs(map, mapExecutor, inputs);

                    // send off the reduce job, and block until we get at least one reduce result
                    sum = reduceSerial(reduce, mapQueue, sum);
                    debugPrint("  Done with cycle of map/reduce");

                    if ( progressFunction != null ) progressFunction.progress(inputs.get(inputs.size()-1));
                } else {
                    // we must be done
                    if ( ! done ) throw new IllegalStateException("Inputs empty but not done");
                }
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            } catch (ExecutionException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            }
        }

        return sum;
    }

    @Requires({"reduce != null", "! mapQueue.isEmpty()"})
    private ReduceType reduceSerial(final NanoSchedulerReduceFunction<MapType, ReduceType> reduce,
                                    final Queue<Future<MapType>> mapQueue,
                                    final ReduceType initSum)
            throws InterruptedException, ExecutionException {
        ReduceType sum = initSum;

        // while mapQueue has something in it to reduce
        for ( final Future<MapType> future : mapQueue ) {
            final MapType value = future.get(); // block until we get the values for this task

            if ( TIME_CALLS ) reduceTimer.restart();
            sum = reduce.apply(value, sum);
            if ( TIME_CALLS ) reduceTimer.stop();
        }

        return sum;
    }

    /**
     * Read up to inputBufferSize elements from inputReader
     *
     * @return a queue of input read in, containing one or more values of InputType read in
     */
    @Requires("inputReader != null")
    @Ensures("result != null")
    private Pair<List<InputType>, Boolean> readInputs(final BlockingQueue<InputDatum> inputReader) throws InterruptedException {
        int n = 0;
        final List<InputType> inputs = new LinkedList<InputType>();
        boolean done = false;

        while ( ! done && n < getBufferSize() ) {
            final InputDatum input = inputReader.take();
            done = input.isLast();
            if ( ! done ) {
                inputs.add(input.datum);
                n++;
            }
        }

        return new Pair<List<InputType>, Boolean>(inputs, done);
    }

    private class InputProducer implements Runnable {
        final Iterator<InputType> inputReader;
        final BlockingQueue<InputDatum> outputQueue;

        public InputProducer(final Iterator<InputType> inputReader, final BlockingQueue<InputDatum> outputQueue) {
            this.inputReader = inputReader;
            this.outputQueue = outputQueue;
        }

        public void run() {
            try {
                while ( inputReader.hasNext() ) {
                    if ( TIME_CALLS ) inputTimer.restart();
                    final InputType input = inputReader.next();
                    if ( TIME_CALLS ) inputTimer.stop();
                    outputQueue.put(new InputDatum(input));
                }

                // add the EOF object so we know we are done
                outputQueue.put(new InputDatum());
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
            }
        }
    }

    private class InputDatum {
        final boolean isLast;
        final InputType datum;

        private InputDatum(final InputType datum) {
            isLast = false;
            this.datum = datum;
        }

        private InputDatum() {
            isLast = true;
            this.datum = null;
        }

        public boolean isLast() {
            return isLast;
        }
    }

    @Requires({"map != null", "! inputs.isEmpty()"})
    private Queue<Future<MapType>> submitMapJobs(final NanoSchedulerMapFunction<InputType, MapType> map,
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
        final NanoSchedulerMapFunction<InputType, MapType> map;

        @Requires({"map != null"})
        private CallableMap(final NanoSchedulerMapFunction<InputType, MapType> map, final InputType inputs) {
            this.input = inputs;
            this.map = map;
        }

        @Override public MapType call() throws Exception {
            if ( TIME_CALLS ) mapTimer.restart();
            if ( debug ) debugPrint("\t\tmap " + input);
            final MapType result = map.apply(input);
            if ( TIME_CALLS ) mapTimer.stop();
            return result;
        }
    }
}
