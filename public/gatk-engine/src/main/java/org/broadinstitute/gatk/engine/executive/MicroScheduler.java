/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.executive;

import com.google.java.contract.Ensures;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.ReadMetrics;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.engine.iterators.NullSAMIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.engine.traversals.*;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.AutoFormattingTime;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.progressmeter.ProgressMeter;
import org.broadinstitute.gatk.utils.threading.ThreadEfficiencyMonitor;

import javax.management.JMException;
import javax.management.MBeanServer;
import javax.management.ObjectName;
import java.io.File;
import java.lang.management.ManagementFactory;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 *
 * General base class for all scheduling algorithms
 * Shards and schedules data in manageable chunks.
 *
 * Creates N TraversalEngines for each data thread for the MicroScheduler.  This is necessary
 * because in the HMS case you have multiple threads executing a traversal engine independently, and
 * these engines may need to create separate resources for efficiency or implementation reasons.  For example,
 * the nanoScheduler creates threads to implement the traversal, and this creation is instance specific.
 * So each HMS thread needs to have it's own distinct copy of the traversal engine if it wants to have
 * N data threads x M nano threads => N * M threads total.  These are borrowed from this microscheduler
 * and returned when done.  Also allows us to tracks all created traversal engines so this microscheduler
 * can properly shut them all down when the scheduling is done.
 *
 */
public abstract class MicroScheduler implements MicroSchedulerMBean {
    protected static final Logger logger = Logger.getLogger(MicroScheduler.class);

    /**
     * The list of all Traversal engines we've created in this micro scheduler
     */
    final List<TraversalEngine> allCreatedTraversalEngines = new LinkedList<TraversalEngine>();

    /**
     * All available engines.  Engines are borrowed and returned when a subclass is actually
     * going to execute the engine on some data.  This allows us to have N copies for
     * N data parallel executions, but without the dangerous code of having local
     * ThreadLocal variables.
     */
    final LinkedList<TraversalEngine> availableTraversalEngines = new LinkedList<TraversalEngine>();

    /**
     * Engines that have been allocated to a key already.
     */
    final HashMap<Object, TraversalEngine> allocatedTraversalEngines = new HashMap<Object, TraversalEngine>();

    /**
     * Counts the number of instances of the class that are currently alive.
     */
    private static int instanceNumber = 0;

    /**
     * The engine invoking this scheduler.
     */
    protected final GenomeAnalysisEngine engine;

    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;
    protected final Collection<ReferenceOrderedDataSource> rods;

    private final MBeanServer mBeanServer;
    private final ObjectName mBeanName;

    /**
     * Threading efficiency monitor for tracking the resource utilization of the GATK
     *
     * may be null
     */
    ThreadEfficiencyMonitor threadEfficiencyMonitor = null;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     *
     * @param walker        Which walker to use.
     * @param reads         the informations associated with the reads
     * @param reference     the reference file
     * @param rods          the rods to include in the traversal
     * @param threadAllocation Number of threads to utilize.
     *
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, ThreadAllocation threadAllocation) {
        if ( threadAllocation.isRunningInParallelMode() ) {
            logger.info(String.format("Running the GATK in parallel mode with %d total threads, " +
                    "%d CPU thread(s) for each of %d data thread(s), of %d processors available on this machine",
                    threadAllocation.getTotalNumThreads(),
                    threadAllocation.getNumCPUThreadsPerDataThread(),
                    threadAllocation.getNumDataThreads(),
                    Runtime.getRuntime().availableProcessors()));
            if ( threadAllocation.getTotalNumThreads() > Runtime.getRuntime().availableProcessors() )
                logger.warn(String.format("Number of requested GATK threads %d is more than the number of " +
                        "available processors on this machine %d", threadAllocation.getTotalNumThreads(),
                        Runtime.getRuntime().availableProcessors()));
        }

        if ( threadAllocation.getNumDataThreads() > 1 ) {
            if (walker.isReduceByInterval())
                throw new UserException.BadArgumentValue("nt", String.format("This run of %s is set up to aggregate results by interval.  Due to a current limitation of the GATK, analyses of this type do not currently support parallel execution.  Please run your analysis without the -nt option or check if this tool has an option to disable per-interval calculations.", engine.getWalkerName(walker.getClass())));

            if ( ! (walker instanceof TreeReducible) ) {
                throw badNT("nt", engine, walker);
            }
        }

        if ( threadAllocation.getNumCPUThreadsPerDataThread() > 1 && ! (walker instanceof NanoSchedulable) ) {
            throw badNT("nct", engine, walker);
        }

        if ( threadAllocation.getNumDataThreads() > 1 ) {
            return new HierarchicalMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
        } else {
            return new LinearMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
        }
    }

    private static UserException badNT(final String parallelArg, final GenomeAnalysisEngine engine, final Walker walker) {
        throw new UserException.BadArgumentValue(parallelArg,
                String.format("The analysis %s currently does not support parallel execution with %s.  " +
                        "Please run your analysis without the %s option.", engine.getWalkerName(walker.getClass()), parallelArg, parallelArg));
    }

    /**
     * Create a microscheduler given the reads and reference.
     *
     * @param walker  the walker to execute with
     * @param reads   The reads.
     * @param reference The reference.
     * @param rods    the rods to include in the traversal
     * @param threadAllocation the allocation of threads to use in the underlying traversal
     */
    protected MicroScheduler(final GenomeAnalysisEngine engine,
                             final Walker walker,
                             final SAMDataSource reads,
                             final IndexedFastaSequenceFile reference,
                             final Collection<ReferenceOrderedDataSource> rods,
                             final ThreadAllocation threadAllocation) {
        this.engine = engine;
        this.reads = reads;
        this.reference = reference;
        this.rods = rods;

        final File progressLogFile = engine.getArguments() == null ? null : engine.getArguments().performanceLog;

        // Creates uninitialized TraversalEngines appropriate for walker and threadAllocation,
        // and adds it to the list of created engines for later shutdown.
        for ( int i = 0; i < threadAllocation.getNumDataThreads(); i++ ) {
            final TraversalEngine traversalEngine = createTraversalEngine(walker, threadAllocation);
            allCreatedTraversalEngines.add(traversalEngine);
            availableTraversalEngines.add(traversalEngine);
        }

        // Create the progress meter, and register it with the analysis engine
        engine.registerProgressMeter(new ProgressMeter(progressLogFile,
                availableTraversalEngines.peek().getTraversalUnits(),
                engine.getRegionsOfGenomeBeingProcessed()));

        // Now that we have a progress meter, go through and initialize the traversal engines
        for ( final TraversalEngine traversalEngine : allCreatedTraversalEngines )
            traversalEngine.initialize(engine, walker, engine.getProgressMeter());

        // JMX does not allow multiple instances with the same ObjectName to be registered with the same platform MXBean.
        // To get around this limitation and since we have no job identifier at this point, register a simple counter that
        // will count the number of instances of this object that have been created in this JVM.
        int thisInstance = instanceNumber++;
        mBeanServer = ManagementFactory.getPlatformMBeanServer();
        try {
            mBeanName = new ObjectName("org.broadinstitute.gatk.engine.executive:type=MicroScheduler,instanceNumber="+thisInstance);
            mBeanServer.registerMBean(this, mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedGATKException("Unable to register microscheduler with JMX", ex);
        }
    }

    /**
     * Really make us a traversal engine of the appropriate type for walker and thread allocation
     *
     * @return a non-null uninitialized traversal engine
     */
    @Ensures("result != null")
    private TraversalEngine createTraversalEngine(final Walker walker, final ThreadAllocation threadAllocation) {
        if (walker instanceof ReadWalker) {
            return new TraverseReadsNano(threadAllocation.getNumCPUThreadsPerDataThread());
        } else if (walker instanceof LocusWalker) {
            return new TraverseLociNano(threadAllocation.getNumCPUThreadsPerDataThread());
        } else if (walker instanceof DuplicateWalker) {
            return new TraverseDuplicates();
        } else if (walker instanceof ReadPairWalker) {
            return new TraverseReadPairs();
        } else if (walker instanceof ActiveRegionWalker) {
            return new TraverseActiveRegions(threadAllocation.getNumCPUThreadsPerDataThread());
        } else {
            throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
        }
    }


    /**
     * Return the ThreadEfficiencyMonitor we are using to track our resource utilization, if there is one
     *
     * @return the monitor, or null if none is active
     */
    public ThreadEfficiencyMonitor getThreadEfficiencyMonitor() {
        return threadEfficiencyMonitor;
    }

    /**
     * Inform this Microscheduler to use the efficiency monitor used to create threads in subclasses
     *
     * @param threadEfficiencyMonitor
     */
    public void setThreadEfficiencyMonitor(final ThreadEfficiencyMonitor threadEfficiencyMonitor) {
        this.threadEfficiencyMonitor = threadEfficiencyMonitor;
    }

    /**
     * Should we stop all execution work and exit gracefully?
     *
     * Returns true in the case where some external signal or time limit has been received, indicating
     * that this GATK shouldn't continue executing.  This isn't a kill signal, it is really a "shutdown
     * gracefully at the next opportunity" signal.  Concrete implementations of the MicroScheduler
     * examine this value as often as reasonable and, if it returns true, stop what they are doing
     * at the next available opportunity, shutdown their resources, call notify done, and return.
     *
     * @return true if we should abort execution, or false otherwise
     */
    protected boolean abortExecution() {
        final boolean abort = engine.exceedsRuntimeLimit();
        if ( abort ) {
            final AutoFormattingTime aft = new AutoFormattingTime(engine.getRuntimeLimitInNanoseconds(), -1, 4);
            logger.info("Aborting execution (cleanly) because the runtime has exceeded the requested maximum " + aft);
        }
        return abort;
    }

    /**
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, Iterable<Shard> shardStrategy);

    /**
     * Tells this MicroScheduler that the execution of one of the subclass of this object as started
     *
     * Must be called when the implementation of execute actually starts up
     *
     * Currently only starts the progress meter timer running, but other start up activities could be incorporated
     */
    protected void startingExecution() {
        engine.getProgressMeter().start();
    }

    /**
     * Retrieves the object responsible for tracking and managing output.
     * @return An output tracker, for loading data in and extracting results.  Will not be null.
     */
    public abstract OutputTracker getOutputTracker();

    /**
     * Gets the an iterator over the given reads, which will iterate over the reads in the given shard.
     * @param shard the shard to use when querying reads.
     * @return an iterator over the reads specified in the shard.
     */
    protected GATKSAMIterator getReadIterator(Shard shard) {
        return (!reads.isEmpty()) ? reads.seek(shard) : new NullSAMIterator();
    }

    /**
     * Must be called by subclasses when execute is done
     */
    protected void executionIsDone() {
        engine.getProgressMeter().notifyDone(engine.getCumulativeMetrics().getNumIterations());
        printReadFilteringStats();
        shutdownTraversalEngines();

        // Print out the threading efficiency of this HMS, if state monitoring is enabled
        if ( threadEfficiencyMonitor != null ) {
            // include the master thread information
            threadEfficiencyMonitor.threadIsDone(Thread.currentThread());
            threadEfficiencyMonitor.printUsageInformation(logger);
        }
    }

    /**
     * Shutdown all of the created engines, and clear the list of created engines, dropping
     * pointers to the traversal engines
     */
    public synchronized void shutdownTraversalEngines() {
        for ( final TraversalEngine te : allCreatedTraversalEngines)
            te.shutdown();

        allCreatedTraversalEngines.clear();
        availableTraversalEngines.clear();
    }

    /**
     * Prints out information about number of reads observed and filtering, if any reads were used in the traversal
     *
     * Looks like:
     *
     * INFO  10:40:47,370 MicroScheduler - 22 reads were filtered out during traversal out of 101 total (21.78%)
     * INFO  10:40:47,370 MicroScheduler -   -> 1 reads (0.99% of total) failing BadMateFilter
     * INFO  10:40:47,370 MicroScheduler -   -> 20 reads (19.80% of total) failing DuplicateReadFilter
     * INFO  10:40:47,370 MicroScheduler -   -> 1 reads (0.99% of total) failing FailsVendorQualityCheckFilter
     */
    private void printReadFilteringStats() {
        final ReadMetrics cumulativeMetrics = engine.getCumulativeMetrics();
        if ( cumulativeMetrics.getNumReadsSeen() > 0 ) {
            // count up the number of skipped reads by summing over all filters
            long nSkippedReads = 0L;
            for ( final long countsByFilter : cumulativeMetrics.getCountsByFilter().values())
                nSkippedReads += countsByFilter;

            logger.info(String.format("%d reads were filtered out during the traversal out of approximately %d total reads (%.2f%%)",
                    nSkippedReads,
                    cumulativeMetrics.getNumReadsSeen(),
                    100.0 * MathUtils.ratio(nSkippedReads, cumulativeMetrics.getNumReadsSeen())));

            for ( final Map.Entry<String, Long> filterCounts : cumulativeMetrics.getCountsByFilter().entrySet() ) {
                long count = filterCounts.getValue();
                logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                        count, 100.0 * MathUtils.ratio(count,cumulativeMetrics.getNumReadsSeen()), filterCounts.getKey()));
            }
        }
    }

    /**
     * Gets the engine that created this microscheduler.
     * @return The engine owning this microscheduler.
     */
    public GenomeAnalysisEngine getEngine() { return engine; }

    /**
     * Returns data source maintained by this scheduler
     * @return
     */
    public SAMDataSource getSAMDataSource() { return reads; }

    /**
     * Returns the reference maintained by this scheduler.
     * @return The reference maintained by this scheduler.
     */
    public IndexedFastaSequenceFile getReference() { return reference; }

    protected void cleanup() {
        try {
            mBeanServer.unregisterMBean(mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedGATKException("Unable to unregister microscheduler with JMX", ex);
        }
    }

    /**
     * Returns a traversal engine suitable for use, associated with key
     *
     * Key is an arbitrary object that is used to retrieve the same traversal
     * engine over and over.  This can be important in the case where the
     * traversal engine has data associated with it in some other context,
     * and we need to ensure that the context always sees the same traversal
     * engine.  This happens in the HierarchicalMicroScheduler, where you want
     * the a thread executing traversals to retrieve the same engine each time,
     * as outputs are tracked w.r.t. that engine.
     *
     * If no engine is associated with key yet, pops the next available engine
     * from the available ones maintained by this
     * microscheduler.  Note that it's a runtime error to pop a traversal engine
     * from this scheduler if there are none available.  Callers that
     * once pop'd an engine for use must return it with returnTraversalEngine
     *
     * @param key the key to associate with this engine
     * @return a non-null TraversalEngine suitable for execution in this scheduler
     */
    @Ensures("result != null")
    protected synchronized TraversalEngine borrowTraversalEngine(final Object key) {
        if ( key == null ) throw new IllegalArgumentException("key cannot be null");

        final TraversalEngine engine = allocatedTraversalEngines.get(key);
        if ( engine == null ) {
            if ( availableTraversalEngines.isEmpty() )
                throw new IllegalStateException("no traversal engines were available");
            allocatedTraversalEngines.put(key, availableTraversalEngines.pop());
            return allocatedTraversalEngines.get(key);
        } else {
            return engine;
        }
    }

    /**
     * Return a borrowed traversal engine to this MicroScheduler, for later use
     * in another traversal execution
     *
     * @param key the key used to id the engine, provided to the borrowTraversalEngine function
     * @param traversalEngine the borrowed traversal engine.  Must have been previously borrowed.
     */
    protected synchronized void returnTraversalEngine(final Object key, final TraversalEngine traversalEngine) {
        if ( traversalEngine == null )
            throw new IllegalArgumentException("Attempting to push a null traversal engine");
        if ( ! allCreatedTraversalEngines.contains(traversalEngine) )
            throw new IllegalArgumentException("Attempting to push a traversal engine not created by this MicroScheduler" + engine);
        if ( ! allocatedTraversalEngines.containsKey(key) )
            throw new IllegalArgumentException("No traversal engine was never checked out with key " + key);

        // note there's nothing to actually do here, but a function implementation
        // might want to do something
    }
}
