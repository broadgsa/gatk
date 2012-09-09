/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.executive;

import com.google.java.contract.Ensures;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.iterators.NullSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.progressmeter.ProgressMeter;
import org.broadinstitute.sting.utils.threading.ThreadEfficiencyMonitor;

import javax.management.JMException;
import javax.management.MBeanServer;
import javax.management.ObjectName;
import java.io.File;
import java.lang.management.ManagementFactory;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 *
 * General base class for all scheduling algorithms
 */

/** Shards and schedules data in manageable chunks. */
public abstract class MicroScheduler implements MicroSchedulerMBean {
    // TODO -- remove me and retire non nano scheduled versions of traversals
    private final static boolean USE_NANOSCHEDULER_FOR_EVERYTHING = true;
    protected static final Logger logger = Logger.getLogger(MicroScheduler.class);

    /**
     * Counts the number of instances of the class that are currently alive.
     */
    private static int instanceNumber = 0;

    /**
     * The engine invoking this scheduler.
     */
    protected final GenomeAnalysisEngine engine;

    private final TraversalEngineCreator traversalEngineCreator;
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

    final ProgressMeter progressMeter;

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
            logger.info(String.format("Running the GATK in parallel mode with %d CPU thread(s) for each of %d data thread(s)",
                    threadAllocation.getNumCPUThreadsPerDataThread(), threadAllocation.getNumDataThreads()));
        }

        if ( threadAllocation.getNumDataThreads() > 1 ) {
            if (walker.isReduceByInterval())
                throw new UserException.BadArgumentValue("nt", String.format("The analysis %s aggregates results by interval.  Due to a current limitation of the GATK, analyses of this type do not currently support parallel execution.  Please run your analysis without the -nt option.", engine.getWalkerName(walker.getClass())));

            if ( ! (walker instanceof TreeReducible) ) {
                throw badNT("nt", engine, walker);
            } else {
                return new HierarchicalMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
            }
        } else {
            if ( threadAllocation.getNumCPUThreadsPerDataThread() > 1 && ! (walker instanceof NanoSchedulable) )
                throw badNT("nct", engine, walker);
            return new LinearMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
        }
    }

    private static UserException badNT(final String parallelArg, final GenomeAnalysisEngine engine, final Walker walker) {
        throw new UserException.BadArgumentValue("nt",
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
        this.traversalEngineCreator = new TraversalEngineCreator(walker, threadAllocation);

        final File progressLogFile = engine.getArguments() == null ? null : engine.getArguments().performanceLog;
        this.progressMeter = new ProgressMeter(progressLogFile,
                traversalEngineCreator.getTraversalUnits(),
                engine.getRegionsOfGenomeBeingProcessed());

        // JMX does not allow multiple instances with the same ObjectName to be registered with the same platform MXBean.
        // To get around this limitation and since we have no job identifier at this point, register a simple counter that
        // will count the number of instances of this object that have been created in this JVM.
        int thisInstance = instanceNumber++;
        mBeanServer = ManagementFactory.getPlatformMBeanServer();
        try {
            mBeanName = new ObjectName("org.broadinstitute.sting.gatk.executive:type=MicroScheduler,instanceNumber="+thisInstance);
            mBeanServer.registerMBean(this, mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedStingException("Unable to register microscheduler with JMX", ex);
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
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, Iterable<Shard> shardStrategy);

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
    protected StingSAMIterator getReadIterator(Shard shard) {
        return (!reads.isEmpty()) ? reads.seek(shard) : new NullSAMIterator();
    }

    /**
     * Must be called by subclasses when execute is done
     */
    protected void executionIsDone() {
        progressMeter.notifyDone(engine.getCumulativeMetrics().getNumIterations());
        printReadFilteringStats();

        for ( final TraversalEngine te : traversalEngineCreator.getCreatedEngines() )
            te.shutdown();

        // Print out the threading efficiency of this HMS, if state monitoring is enabled
        if ( threadEfficiencyMonitor != null ) {
            // include the master thread information
            threadEfficiencyMonitor.threadIsDone(Thread.currentThread());
            threadEfficiencyMonitor.printUsageInformation(logger);
        }
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

            logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
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
            throw new ReviewedStingException("Unable to unregister microscheduler with JMX", ex);
        }
    }

    /**
     * Returns a traversal engine suitable for use in this thread.
     *
     * May create a new traversal engine for this thread, if this is the first
     * time this thread ever asked for a TraversalEngine.
     *
     * @return a non-null TraversalEngine suitable for execution in this scheduler
     */
    public TraversalEngine getTraversalEngine() {
        return traversalEngineCreator.get();
    }

    /**
     * ThreadLocal TraversalEngine creator
     *
     * TraversalEngines are thread local variables to the MicroScheduler.  This is necessary
     * because in the HMS case you have multiple threads executing a traversal engine independently, and
     * these engines may need to create separate resources for efficiency or implementation reasons.  For example,
     * the nanoScheduler creates threads to implement the traversal, and this creation is instance specific.
     * So each HMS thread needs to have it's own distinct copy of the traversal engine if it wants to have
     * N data threads x M nano threads => N * M threads total.
     *
     * This class also tracks all created traversal engines so this microscheduler can properly
     * shut them all down when the scheduling is done.
     */
    private class TraversalEngineCreator extends ThreadLocal<TraversalEngine> {
        final List<TraversalEngine> createdEngines = new LinkedList<TraversalEngine>();
        final Walker walker;
        final ThreadAllocation threadAllocation;

        /**
         * Creates an initialized TraversalEngine appropriate for walker and threadAllocation,
         * and adds it to the list of created engines for later shutdown.
         *
         * @return a non-null traversal engine
         */
        @Override
        protected synchronized TraversalEngine initialValue() {
            final TraversalEngine traversalEngine = createEngine();
            traversalEngine.initialize(engine, progressMeter);
            createdEngines.add(traversalEngine);
            return traversalEngine;
        }

        /**
         * Returns the traversal units for traversal engines created here.
         *
         * This (unfortunately) creates an uninitialized tmp. TraversalEngine so we can get
         * it's traversal units, and then immediately shuts it down...
         *
         * @return the traversal unit as returned by getTraversalUnits of TraversalEngines created here
         */
        protected String getTraversalUnits() {
            final TraversalEngine tmp = createEngine();
            final String units = tmp.getTraversalUnits();
            tmp.shutdown();
            return units;
        }

        /**
         * Really make us a traversal engine of the appropriate type for walker and thread allocation
         *
         * @return a non-null uninitialized traversal engine
         */
        @Ensures("result != null")
        protected TraversalEngine createEngine() {
            if (walker instanceof ReadWalker) {
                if ( USE_NANOSCHEDULER_FOR_EVERYTHING || threadAllocation.getNumCPUThreadsPerDataThread() > 1 )
                    return new TraverseReadsNano(threadAllocation.getNumCPUThreadsPerDataThread());
                else
                    return new TraverseReads();
            } else if (walker instanceof LocusWalker) {
                if ( USE_NANOSCHEDULER_FOR_EVERYTHING || threadAllocation.getNumCPUThreadsPerDataThread() > 1 )
                    return new TraverseLociNano(threadAllocation.getNumCPUThreadsPerDataThread());
                else
                    return new TraverseLociLinear();
            } else if (walker instanceof DuplicateWalker) {
                return new TraverseDuplicates();
            } else if (walker instanceof ReadPairWalker) {
                return new TraverseReadPairs();
            } else if (walker instanceof ActiveRegionWalker) {
                return new TraverseActiveRegions();
            } else {
                throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
            }
        }

        /**
         * Create a TraversalEngineCreator that makes TraversalEngines appropriate for walker and threadAllocation
         *
         * @param walker the walker we need traversal engines for
         * @param threadAllocation what kind of threading will we use in the traversal?
         */
        @com.google.java.contract.Requires({"walker != null", "threadAllocation != null"})
        public TraversalEngineCreator(final Walker walker, final ThreadAllocation threadAllocation) {
            super();
            this.walker = walker;
            this.threadAllocation = threadAllocation;
        }

        /**
         * Get the list of all traversal engines we've created
         *
         * @return a non-null list of engines created so far
         */
        @Ensures("result != null")
        public List<TraversalEngine> getCreatedEngines() {
            return createdEngines;
        }
    }
}
