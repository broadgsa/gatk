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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.NullSAMIterator;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.ReadMetrics;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.util.*;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.threading.GenomeLocProcessingTracker;
import org.broadinstitute.sting.utils.threading.NoOpGenomeLocProcessingTracker;
import org.broadinstitute.sting.utils.threading.SharedFileGenomeLocProcessingTracker;

import javax.management.JMException;
import javax.management.MBeanServer;
import javax.management.ObjectName;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 * To change this template use File | Settings | File Templates.
 */

/** Shards and schedules data in manageable chunks. */
public abstract class MicroScheduler implements MicroSchedulerMBean {
    protected static Logger logger = Logger.getLogger(MicroScheduler.class);

    /**
     * Counts the number of instances of the class that are currently alive.
     */
    private static int instanceNumber = 0;

    /**
     * The engine invoking this scheduler.
     */
    protected final GenomeAnalysisEngine engine;

    protected final TraversalEngine traversalEngine;
    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;
    protected final Collection<ReferenceOrderedDataSource> rods;

    private final MBeanServer mBeanServer;
    private final ObjectName mBeanName;

    private GenomeLocProcessingTracker processingTracker;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     *
     * @param walker        Which walker to use.
     * @param reads         the informations associated with the reads
     * @param reference     the reference file
     * @param rods          the rods to include in the traversal
     * @param nThreadsToUse Number of threads to utilize.
     *
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, int nThreadsToUse) {
        if (walker instanceof TreeReducible && nThreadsToUse > 1) {
            if(walker.isReduceByInterval())
                throw new UserException.BadArgumentValue("nt", String.format("The analysis %s aggregates results by interval.  Due to a current limitation of the GATK, analyses of this type do not currently support parallel execution.  Please run your analysis without the -nt option.", engine.getWalkerName(walker.getClass())));
            if(walker instanceof ReadWalker)
                throw new UserException.BadArgumentValue("nt", String.format("The analysis %s is a read walker.  Due to a current limitation of the GATK, analyses of this type do not currently support parallel execution.  Please run your analysis without the -nt option.", engine.getWalkerName(walker.getClass())));
            logger.info(String.format("Running the GATK in parallel mode with %d concurrent threads",nThreadsToUse));
            return new HierarchicalMicroScheduler(engine, walker, reads, reference, rods, nThreadsToUse);
        } else {
            if(nThreadsToUse > 1)
                throw new UserException.BadArgumentValue("nt", String.format("The analysis %s currently does not support parallel execution.  Please run your analysis without the -nt option.", engine.getWalkerName(walker.getClass())));
            return new LinearMicroScheduler(engine, walker, reads, reference, rods);
        }
    }

    /**
     * Create a microscheduler given the reads and reference.
     *
     * @param walker  the walker to execute with
     * @param reads   The reads.
     * @param reference The reference.
     * @param rods    the rods to include in the traversal
     */
    protected MicroScheduler(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        this.engine = engine;
        this.reads = reads;
        this.reference = reference;
        this.rods = rods;

        if (walker instanceof ReadWalker) {
            traversalEngine = new TraverseReads();
        } else if (walker instanceof LocusWalker) {
            traversalEngine = new TraverseLoci();
        } else if (walker instanceof DuplicateWalker) {
            traversalEngine = new TraverseDuplicates();
        } else if (walker instanceof ReadPairWalker) {
            traversalEngine = new TraverseReadPairs();
        } else {
            throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
        }        

        traversalEngine.initialize(engine);

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

        // create the processing tracker
        if ( engine.getArguments().processingTrackerFile != null ) {
            if ( engine.getArguments().restartProcessingTracker && engine.getArguments().processingTrackerFile.exists() ) {
                engine.getArguments().processingTrackerFile.delete();
                logger.info("Deleting ProcessingTracker file " + engine.getArguments().processingTrackerFile);
            }

            processingTracker = new SharedFileGenomeLocProcessingTracker(engine.getArguments().processingTrackerFile, engine.getGenomeLocParser());
            logger.info("Creating ProcessingTracker using shared file " + engine.getArguments().processingTrackerFile);
        } else {
            processingTracker = new NoOpGenomeLocProcessingTracker();
        }
    }

    /**
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, ShardStrategy shardStrategy);

    protected boolean claimShard(Shard shard) {
        if ( shard.getGenomeLocs() == null ) {
            if ( engine.getArguments().processingTrackerFile != null )
                throw new UserException.BadArgumentValue("processingTrackerFile", "Cannot use processing tracking with unindexed data");
            return true;
        } else {
            GenomeLoc shardSpan = shardSpan(shard);

            GenomeLocProcessingTracker.ProcessingLoc proc = processingTracker.claimOwnership(shardSpan, engine.getName());
            boolean actuallyProcess = proc.isOwnedBy(engine.getName());
            //logger.debug(String.format("Shard %s claimed by %s => owned by me %b", shard, proc.getOwner(), actuallyProcess));

            if ( ! actuallyProcess )
                logger.info(String.format("DISTRIBUTED GATK: Shard %s already processed by %s", shard, proc.getOwner()));

            return actuallyProcess;
        }
    }

    private GenomeLoc shardSpan(Shard shard) {
        if ( shard == null ) throw new ReviewedStingException("Shard is null!");
        int start = Integer.MAX_VALUE;
        int stop = Integer.MIN_VALUE;
        String contig = null;

        for ( GenomeLoc loc : shard.getGenomeLocs() ) {
            if ( GenomeLoc.isUnmapped(loc) )
                // special case the unmapped region marker, just abort out
                return loc;
            contig = loc.getContig();
            if ( loc.getStart() < start ) start = loc.getStart();
            if ( loc.getStop() > stop ) stop = loc.getStop();
        }
        return engine.getGenomeLocParser().createGenomeLoc(contig, start, stop);
    }

    // todo -- the execution code in the schedulers is duplicated and slightly different -- should be merged
//    protected boolean executeShard(Walker walker, Shard shard, Accumulator accumulator) {
//        if(shard.getShardType() == Shard.ShardType.LOCUS) {
//            LocusWalker lWalker = (LocusWalker)walker;
//
//            WindowMaker windowMaker = new WindowMaker(shard, engine.getGenomeLocParser(), getReadIterator(shard), shard.getGenomeLocs(), lWalker.getDiscards(), engine.getSampleMetadata());
//
//            // ShardTraverser
//
//            ShardDataProvider dataProvider = null;
//            for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
//                dataProvider = new LocusShardDataProvider(shard,iterator.getSourceInfo(),engine.getGenomeLocParser(),iterator.getLocus(),iterator,reference,rods);
//                Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
//                accumulator.accumulate(dataProvider,result);
//                dataProvider.close();
//            }
//            if (dataProvider != null) dataProvider.close();
//            windowMaker.close();
//        }
//        else {
//            ShardDataProvider dataProvider = new ReadShardDataProvider(shard,engine.getGenomeLocParser(),getReadIterator(shard),reference,rods);
//            Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
//            accumulator.accumulate(dataProvider,result);
//            dataProvider.close();
//        }
//
//
//        // ShardTraverser
//
//        for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
//            accumulator = traversalEngine.traverse( walker, dataProvider, accumulator );
//            dataProvider.close();
//        }
//
//        return true;
//    }

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
     * Print summary information for the analysis.
     * @param sum The final reduce output.
     */
    protected void printOnTraversalDone(Object sum, ReadMetrics metrics) {
        traversalEngine.printOnTraversalDone(metrics);
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

    /**
     * Gets the filename to which performance data is currently being written.
     * @return Filename to which performance data is currently being written.
     */
    public String getPerformanceLogFileName() {
        return traversalEngine.getPerformanceLogFileName();
    }

    /**
     * Set the filename of the log for performance.  If set,
     * @param fileName filename to use when writing performance data.
     */
    public void setPerformanceLogFileName(String fileName) {
        traversalEngine.setPerformanceLogFileName(fileName);
    }

    /**
     * Gets the frequency with which performance data is written.
     * @return Frequency, in seconds, of performance log writes.
     */
    public long getPerformanceProgressPrintFrequencySeconds() {
        return traversalEngine.getPerformanceProgressPrintFrequencySeconds();
    }    

    /**
     * How often should the performance log message be written?
     * @param seconds number of seconds between messages indicating performance frequency.
     */
    public void setPerformanceProgressPrintFrequencySeconds(long seconds) {
        traversalEngine.setPerformanceProgressPrintFrequencySeconds(seconds);
    }

    protected void cleanup() {
        try {
            mBeanServer.unregisterMBean(mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedStingException("Unable to unregister microscheduler with JMX", ex);
        }
    }
}
