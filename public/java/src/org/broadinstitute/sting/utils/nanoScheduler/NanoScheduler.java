package org.broadinstitute.sting.utils.nanoScheduler;

import org.apache.log4j.Logger;
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
 * nano <- new Nanoschedule(bufferSize, numberOfMapElementsToProcessTogether, nThreads)
 * List[Input] outerData : outerDataLoop )
 *   result = nano.execute(outerData.iterator(), map, reduce)
 *
 * bufferSize determines how many elements from the input stream are read in one go by the
 * nanoscheduler.  The scheduler may hold up to bufferSize in memory at one time, as well
 * as up to bufferSize map results as well.
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

    private final static int MAP_BUFFER_SIZE_SCALE_FACTOR = 100;

    final int bufferSize;
    final int nThreads;
    final ExecutorService inputExecutor;
    final ExecutorService reduceExecutor;
    final ExecutorService mapExecutor;
    final Semaphore mapQueueSizeManagingSemaphone;

    boolean shutdown = false;
    boolean debug = false;
    private NSProgressFunction<InputType> progressFunction = null;

    /**
     * Tracks the combined runtime profiles across all created nano schedulers
     */
    final static private NSRuntimeProfile combinedNSRuntimeProfiler = new NSRuntimeProfile();

    /**
     * The profile specific to this nano scheduler
     */
    final private NSRuntimeProfile myNSRuntimeProfile = new NSRuntimeProfile();

    /**
     * Create a new nanoscheduler with the desire characteristics requested by the argument
     *
     * @param nThreads the number of threads to use to get work done, in addition to the
     *                 thread calling execute
     */
    public NanoScheduler(final int nThreads) {
        this(nThreads*100, nThreads);
    }

    protected NanoScheduler(final int bufferSize, final int nThreads) {
        if ( bufferSize < 1 ) throw new IllegalArgumentException("bufferSize must be >= 1, got " + bufferSize);
        if ( nThreads < 1 ) throw new IllegalArgumentException("nThreads must be >= 1, got " + nThreads);

        this.bufferSize = bufferSize;
        this.nThreads = nThreads;

        if ( nThreads == 1 ) {
            this.mapExecutor = this.inputExecutor = this.reduceExecutor = null;
            mapQueueSizeManagingSemaphone = null;
        } else {
            this.mapExecutor = Executors.newFixedThreadPool(nThreads, new NamedThreadFactory("NS-map-thread-%d"));
            mapQueueSizeManagingSemaphone = new Semaphore(this.bufferSize);

            this.inputExecutor = Executors.newSingleThreadExecutor(new NamedThreadFactory("NS-input-thread-%d"));
            this.reduceExecutor = Executors.newSingleThreadExecutor(new NamedThreadFactory("NS-reduce-thread-%d"));
        }

        // start timing the time spent outside of the nanoScheduler
        myNSRuntimeProfile.outsideSchedulerTimer.start();
    }

    /**
     * The number of parallel map threads in use with this NanoScheduler
     * @return
     */
//    @Ensures("result > 0")
    public int getnThreads() {
        return nThreads;
    }

    /**
     * The input buffer size used by this NanoScheduler
     * @return
     */
//    @Ensures("result > 0")
    public int getBufferSize() {
        return this.bufferSize;
    }

    /**
     * Tells this nanoScheduler to shutdown immediately, releasing all its resources.
     *
     * After this call, execute cannot be invoked without throwing an error
     */
    public void shutdown() {
        myNSRuntimeProfile.outsideSchedulerTimer.stop();

        // add my timing information to the combined NS runtime profile
        combinedNSRuntimeProfiler.combine(myNSRuntimeProfile);

        if ( nThreads > 1 ) {
            shutdownExecutor("inputExecutor", inputExecutor);
            shutdownExecutor("mapExecutor", mapExecutor);
            shutdownExecutor("reduceExecutor", reduceExecutor);
        }

        shutdown = true;
    }

    public void printRuntimeProfile() {
        myNSRuntimeProfile.log(logger);
    }

    public static void printCombinedRuntimeProfile() {
        if ( combinedNSRuntimeProfiler.totalRuntimeInSeconds() > 0.1 )
            combinedNSRuntimeProfiler.log(logger);
    }

    protected double getTotalRuntime() {
        return myNSRuntimeProfile.totalRuntimeInSeconds();
    }

    /**
     * Helper function to cleanly shutdown an execution service, checking that the execution
     * state is clean when it's done.
     *
     * @param name a string name for error messages for the executorService we are shutting down
     * @param executorService the executorService to shut down
     */
//    @Requires({"name != null", "executorService != null"})
//    @Ensures("executorService.isShutdown()")
    private void shutdownExecutor(final String name, final ExecutorService executorService) {
        if ( executorService.isShutdown() || executorService.isTerminated() )
            throw new IllegalStateException("Executor service " + name + " is already shut down!");

        final List<Runnable> remaining = executorService.shutdownNow();
        if ( ! remaining.isEmpty() )
            throw new IllegalStateException(remaining.size() + " remaining tasks found in an executor " + name + ", unexpected behavior!");
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
//    @Requires("format != null")
    private void debugPrint(final String format, Object ... args) {
        if ( isDebug() )
            logger.warn("Thread " + Thread.currentThread().getId() + ":" + String.format(format, args));
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

        myNSRuntimeProfile.outsideSchedulerTimer.stop();

        ReduceType result;
        if ( ALLOW_SINGLE_THREAD_FASTPATH && getnThreads() == 1 ) {
            result = executeSingleThreaded(inputReader, map, initialValue, reduce);
        } else {
            result = executeMultiThreaded(inputReader, map, initialValue, reduce);
        }

        myNSRuntimeProfile.outsideSchedulerTimer.restart();
        return result;
    }

    /**
     * Simple efficient reference implementation for single threaded execution.
     *
     * @return the reduce result of this map/reduce job
     */
//    @Requires({"inputReader != null", "map != null", "reduce != null"})
    private ReduceType executeSingleThreaded(final Iterator<InputType> inputReader,
                                             final NSMapFunction<InputType, MapType> map,
                                             final ReduceType initialValue,
                                             final NSReduceFunction<MapType, ReduceType> reduce) {
        ReduceType sum = initialValue;
        int i = 0;

        while ( true ) {
            // start timer to ensure that both hasNext and next are caught by the timer
            myNSRuntimeProfile.inputTimer.restart();
            if ( ! inputReader.hasNext() ) {
                myNSRuntimeProfile.inputTimer.stop();
                break;
            } else {
                final InputType input = inputReader.next();
                myNSRuntimeProfile.inputTimer.stop();

                // map
                myNSRuntimeProfile.mapTimer.restart();
                final long preMapTime = LOG_MAP_TIMES ? 0 : myNSRuntimeProfile.mapTimer.currentTimeNano();
                final MapType mapValue = map.apply(input);
                if ( LOG_MAP_TIMES ) logger.info("MAP TIME " + (myNSRuntimeProfile.mapTimer.currentTimeNano() - preMapTime));
                myNSRuntimeProfile.mapTimer.stop();

                if ( i++ % this.bufferSize == 0 && progressFunction != null )
                    progressFunction.progress(input);

                // reduce
                myNSRuntimeProfile.reduceTimer.restart();
                sum = reduce.apply(mapValue, sum);
                myNSRuntimeProfile.reduceTimer.stop();
            }
        }

        return sum;
    }

    /**
     * Efficient parallel version of Map/Reduce
     *
     * @return the reduce result of this map/reduce job
     */
//    @Requires({"inputReader != null", "map != null", "reduce != null"})
    private ReduceType executeMultiThreaded(final Iterator<InputType> inputReader,
                                            final NSMapFunction<InputType, MapType> map,
                                            final ReduceType initialValue,
                                            final NSReduceFunction<MapType, ReduceType> reduce) {
//        debugPrint("Executing nanoScheduler");
//
//        // a blocking queue that limits the number of input datum to the requested buffer size
//        final BlockingQueue<InputProducer<InputType>.InputValue> inputQueue
//                = new LinkedBlockingDeque<InputProducer<InputType>.InputValue>(bufferSize);
//
//        // a priority queue that stores up to bufferSize elements
//        // produced by completed map jobs.
//        final BlockingQueue<Future<MapResult<MapType>>> mapResultQueue =
//                new LinkedBlockingDeque<Future<MapResult<MapType>>>(bufferSize);
//
//        // Start running the input reader thread
//        inputExecutor.submit(new InputProducer<InputType>(inputReader, myNSRuntimeProfile.inputTimer, inputQueue));
//
//        // Start running the reducer thread
//        final ReducerThread<MapType, ReduceType> reducer
//                = new ReducerThread<MapType, ReduceType>(reduce, myNSRuntimeProfile.reduceTimer, initialValue, mapResultQueue);
//        final Future<ReduceType> reduceResult = reduceExecutor.submit(reducer);
//
//        try {
//            int numJobs = 0;
//
//            while ( true ) {
//                // block on input
//                final InputProducer<InputType>.InputValue inputEnqueueWrapped = inputQueue.take();
//
//                if ( ! inputEnqueueWrapped.isLast() ) {
//                    // get the object itself
//                    final InputType input = inputEnqueueWrapped.getValue();
//
//                    // the next map call has jobID + 1
//                    numJobs++;
//
//                    // send job for map via the completion service
//                    final CallableMap doMap = new CallableMap(map, numJobs, input);
//                    final Future<MapResult<MapType>> mapJob = mapExecutor.submit(doMap);
//                    mapResultQueue.put(mapJob);
//
//                    debugPrint("  Done with cycle of map/reduce");
//
//                    if ( numJobs % bufferSize == 0 && progressFunction != null )
//                        progressFunction.progress(input);
//                } else {
//                    mapResultQueue.put(new FutureValue<MapResult<MapType>>(new MapResult<MapType>()));
//                    return reduceResult.get(); // wait for our result of reduce
//                }
//            }
//        } catch (InterruptedException ex) {
//            throw new ReviewedStingException("got execution exception", ex);
//        } catch (ExecutionException ex) {
//            throw new ReviewedStingException("got execution exception", ex);
//        }
//    }

        debugPrint("Executing nanoScheduler");

        final InputProducer<InputType> inputProducer =
                new InputProducer<InputType>(inputReader, myNSRuntimeProfile.inputTimer);

        // a priority queue that stores up to bufferSize elements
        // produced by completed map jobs.
        final PriorityBlockingQueue<MapResult<MapType>> mapResultQueue =
                new PriorityBlockingQueue<MapResult<MapType>>();

        final Reducer<MapType, ReduceType> reducer
                = new Reducer<MapType, ReduceType>(reduce, myNSRuntimeProfile.reduceTimer, initialValue);

        try {
            int jobID = -1;

            while ( inputProducer.hasNextNow() ) {
                mapQueueSizeManagingSemaphone.acquire();
                jobID++;
                debugPrint("Submitting job with id %d", jobID);
                mapExecutor.submit(new ReadMapReduceJob(jobID, inputProducer, mapResultQueue, map, reducer));
            }

            debugPrint("Setting last job id %d", jobID);
            reducer.setLastJobID(jobID); // the last actually submitted job id is jobID - 1

            return reducer.waitForFinalReduce();
        } catch (InterruptedException ex) {
            throw new ReviewedStingException("got execution exception", ex);
//        } catch (ExecutionException ex) {
//            throw new ReviewedStingException("got execution exception", ex);
        }
    }

    private class ReadMapReduceJob implements Runnable {
        final int jobID;
        final InputProducer<InputType> inputProducer;
        final BlockingQueue<MapResult<MapType>> mapResultQueue;
        final NSMapFunction<InputType, MapType> map;
        final Reducer<MapType, ReduceType> reducer;

        private ReadMapReduceJob(final int jobID,
                                 final InputProducer<InputType> inputProducer,
                                 final BlockingQueue<MapResult<MapType>> mapResultQueue,
                                 final NSMapFunction<InputType, MapType> map,
                                 final Reducer<MapType, ReduceType> reducer) {
            this.jobID = jobID;
            this.inputProducer = inputProducer;
            this.mapResultQueue = mapResultQueue;
            this.map = map;
            this.reducer = reducer;
        }

        @Override
        public void run() {
            try {
                debugPrint("Running ReadMapReduceJob " + jobID);
                final InputProducer<InputType>.InputValue inputWrapper = inputProducer.next();

                final MapResult<MapType> result;
                if ( ! inputWrapper.isLast() ) {
                    // just skip doing anything if we don't have work to do, which is possible
                    // because we don't necessarily know how much input there is when we queue
                    // up our jobs
                    final InputType input = inputWrapper.getValue();

                    // map
                    myNSRuntimeProfile.mapTimer.restart();
                    final long preMapTime = LOG_MAP_TIMES ? 0 : myNSRuntimeProfile.mapTimer.currentTimeNano();
                    final MapType mapValue = map.apply(input);
                    if ( LOG_MAP_TIMES ) logger.info("MAP TIME " + (myNSRuntimeProfile.mapTimer.currentTimeNano() - preMapTime));
                    myNSRuntimeProfile.mapTimer.stop();

                    // enqueue the result into the mapResultQueue
                    result = new MapResult<MapType>(mapValue, jobID);

                    if ( jobID % bufferSize == 0 && progressFunction != null )
                        progressFunction.progress(input);
                } else {
                    // if there's no input we push empty MapResults with jobIDs for synchronization with Reducer
                    result = new MapResult<MapType>(jobID);
                }

                mapResultQueue.put(result);
                debugPrint("  Pushed MapResult with job id %d", jobID);

                final int nReduced = reducer.reduceAsMuchAsPossible(mapResultQueue);
                debugPrint("  reduced %d values", nReduced);

                // we finished a map job, release the job queue semaphore
                mapQueueSizeManagingSemaphone.release();
            } catch (InterruptedException ex) {
                throw new ReviewedStingException("got execution exception", ex);
//            } catch (ExecutionException ex) {
//                throw new ReviewedStingException("got execution exception", ex);
            }
        }
    }

//    /**
//     * A simple callable version of the map function for use with the executor pool
//     */
//    private class CallableMap implements Callable<MapResult<MapType>> {
//        final int id;
//        final InputType input;
//        final NSMapFunction<InputType, MapType> map;
//
//        @Requires({"map != null"})
//        private CallableMap(final NSMapFunction<InputType, MapType> map,
//                            final int id,
//                            final InputType input) {
//            this.id = id;
//            this.input = input;
//            this.map = map;
//        }
//
//        @Override
//        public MapResult<MapType> call() {
//            if ( debug ) debugPrint("\t\tmap " + input);
//            myNSRuntimeProfile.mapTimer.restart();
//            final MapType result = map.apply(input);
//            myNSRuntimeProfile.mapTimer.stop();
//            return new MapResult<MapType>(result, id);
//        }
//    }
}
