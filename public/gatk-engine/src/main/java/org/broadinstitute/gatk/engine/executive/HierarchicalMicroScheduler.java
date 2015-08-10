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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.TribbleException;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.engine.io.ThreadGroupOutputTracker;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.MultiThreadedErrorTracker;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.threading.ThreadPoolMonitor;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.concurrent.*;

/**
 * A microscheduler that schedules shards according to a tree-like structure.
 * Requires a special walker tagged with a 'TreeReducible' interface.
 */
public class HierarchicalMicroScheduler extends MicroScheduler implements HierarchicalMicroSchedulerMBean, ReduceTree.TreeReduceNotifier {
    /**
     * How many outstanding output merges are allowed before the scheduler stops
     * allowing new processes and starts merging flat-out.
     */
    private static final int MAX_OUTSTANDING_OUTPUT_MERGES = 50;

    /** Manage currently running threads. */
    private ExecutorService threadPool;

    /**
     * A thread local output tracker for managing output per-thread.
     */
    private ThreadGroupOutputTracker outputTracker = new ThreadGroupOutputTracker();

    private final Queue<TreeReduceTask> reduceTasks = new LinkedList<TreeReduceTask>();

    /**
     * An exception that's occurred in this traversal.  If null, no exception has occurred.
     */
    final MultiThreadedErrorTracker errorTracker = new MultiThreadedErrorTracker();

    /**
     * Queue of incoming shards.
     */
    private Iterator<Shard> traversalTasks;

    /**
     * Keep a queue of shard traversals, and constantly monitor it to see what output
     * merge tasks remain.
     * TODO: Integrate this into the reduce tree.
     */
    private final Queue<ShardTraverser> outputMergeTasks = new LinkedList<ShardTraverser>();

    /** How many shard traversals have run to date? */
    private int totalCompletedTraversals = 0;

    /** What is the total time spent traversing shards? */
    private long totalShardTraverseTime = 0;

    /** What is the total time spent tree reducing shard output? */
    private long totalTreeReduceTime = 0;

    /** How many tree reduces have been completed? */
    private long totalCompletedTreeReduces = 0;

    /** What is the total time spent merging output? */
    private long totalOutputMergeTime = 0;

    /**
     * Create a new hierarchical microscheduler to process the given reads and reference.
     *
     * @param walker           the walker used to process the dataset.
     * @param reads            Reads file(s) to process.
     * @param reference        Reference for driving the traversal.
     * @param threadAllocation How should we apply multi-threaded execution?
     */
    protected HierarchicalMicroScheduler(final GenomeAnalysisEngine engine,
                                         final Walker walker,
                                         final SAMDataSource reads,
                                         final IndexedFastaSequenceFile reference,
                                         final Collection<ReferenceOrderedDataSource> rods,
                                         final ThreadAllocation threadAllocation) {
        super(engine, walker, reads, reference, rods, threadAllocation);

        final int nThreadsToUse = threadAllocation.getNumDataThreads();
        if ( threadAllocation.monitorThreadEfficiency() ) {
            throw new UserException.BadArgumentValue("nt", "Cannot monitor thread efficiency with -nt, sorry");
        }

        this.threadPool = Executors.newFixedThreadPool(nThreadsToUse, new UniqueThreadGroupThreadFactory());
    }

    /**
     * Creates threads for HMS each with a unique thread group.  Critical to
     * track outputs via the ThreadGroupOutputTracker.
     */
    private static class UniqueThreadGroupThreadFactory implements ThreadFactory {
        int counter = 0;

        @Override
        public Thread newThread(Runnable r) {
            final ThreadGroup group = new ThreadGroup("HMS-group-" + counter++);
            return new Thread(group, r);
        }
    }

    public Object execute( Walker walker, Iterable<Shard> shardStrategy ) {
        super.startingExecution();

        // Fast fail for walkers not supporting TreeReducible interface.
        if (!( walker instanceof TreeReducible ))
            throw new IllegalArgumentException("The GATK can currently run in parallel only with TreeReducible walkers");

        this.traversalTasks = shardStrategy.iterator();

        final ReduceTree reduceTree = new ReduceTree(this);
        initializeWalker(walker);

        while (! abortExecution() && (isShardTraversePending() || isTreeReducePending())) {
            // Check for errors during execution.
            errorTracker.throwErrorIfPending();

            // Too many files sitting around taking up space?  Merge them.
            if (isMergeLimitExceeded())
                mergeExistingOutput(false);

            // Wait for the next slot in the queue to become free.
            waitForFreeQueueSlot();

            // Pick the next most appropriate task and run it.  In the interest of
            // memory conservation, hierarchical reduces always run before traversals.
            if (isTreeReduceReady())
                queueNextTreeReduce(walker);
            else if (isShardTraversePending())
                queueNextShardTraverse(walker, reduceTree);
        }

        errorTracker.throwErrorIfPending();

        threadPool.shutdown();

        // Merge any lingering output files.  If these files aren't ready,
        // sit around and wait for them, then merge them.
        mergeExistingOutput(true);

        Object result = null;
        try {
            result = reduceTree.getResult().get();
            notifyTraversalDone(walker,result);
        } catch (ReviewedGATKException ex) {
            throw ex;
        } catch ( ExecutionException ex ) {
            // the thread died and we are failing to get the result, rethrow it as a runtime exception
            throw notifyOfTraversalError(ex.getCause());
        } catch (Exception ex) {
            throw new ReviewedGATKException("Unable to retrieve result", ex);
        }

        // do final cleanup operations
        outputTracker.close();
        cleanup();
        executionIsDone();

        return result;
    }

    /**
     * Run the initialize method of the walker.  Ensure that any calls
     * to the output stream will bypass thread local storage and write
     * directly to the output file.
     * @param walker Walker to initialize.
     */
    protected void initializeWalker(Walker walker) {
        outputTracker.bypassThreadLocalStorage(true);
        try {
            walker.initialize();
        }
        finally {
            outputTracker.bypassThreadLocalStorage(false);
        }
    }

    /**
     * Run the initialize method of the walker.  Ensure that any calls
     * to the output stream will bypass thread local storage and write
     * directly to the output file.
     * @param walker Walker to initialize.
     */
    protected void notifyTraversalDone(Walker walker, Object result) {
        outputTracker.bypassThreadLocalStorage(true);
        try {
            walker.onTraversalDone(result);
        }
        finally {
            outputTracker.bypassThreadLocalStorage(false);
        }
    }

    /**
     * @{inheritDoc}
     */
    public OutputTracker getOutputTracker() {
        return outputTracker;
    }

    /**
     * Returns true if there are unscheduled shard traversal waiting to run.
     *
     * @return true if a shard traversal is waiting; false otherwise.
     */
    protected boolean isShardTraversePending() {
        return traversalTasks.hasNext();
    }

    /**
     * Returns true if there are tree reduces that can be run without
     * blocking.
     *
     * @return true if a tree reduce is ready; false otherwise.
     */
    protected boolean isTreeReduceReady() {
        if (reduceTasks.size() == 0)
            return false;
        return reduceTasks.peek().isReadyForReduce();
    }

    /**
     * Returns true if there are tree reduces that need to be run before
     * the computation is complete.  Returns true if any entries are in the queue,
     * blocked or otherwise.
     *
     * @return true if a tree reduce is pending; false otherwise.
     */
    protected boolean isTreeReducePending() {
        return reduceTasks.size() > 0;
    }

    /**
     * Returns whether the maximum number of files is sitting in the temp directory
     * waiting to be merged back in.
     *
     * @return True if the merging needs to take priority.  False otherwise.
     */
    protected boolean isMergeLimitExceeded() {
        int pendingTasks = 0;
        for( ShardTraverser shardTraverse: outputMergeTasks ) {
            if( !shardTraverse.isComplete() )
                break;
            pendingTasks++;
        }
        return (outputMergeTasks.size() >= MAX_OUTSTANDING_OUTPUT_MERGES);
    }

    /**
     * Merging all output that's sitting ready in the OutputMerger queue into
     * the final data streams.
     */
    protected void mergeExistingOutput( boolean wait ) {
        long startTime = System.currentTimeMillis();

//        logger.warn("MergingExistingOutput");
//        printOutputMergeTasks();

        // Create a list of the merge tasks that will be performed in this run of the mergeExistingOutput().
        Queue<ShardTraverser> mergeTasksInSession = new LinkedList<ShardTraverser>();
        while( !outputMergeTasks.isEmpty() ) {
            ShardTraverser traverser = outputMergeTasks.peek();

            // If the next traversal isn't done and we're not supposed to wait, we've found our working set.  Continue.
            if( !traverser.isComplete() && !wait )
                break;

            outputMergeTasks.remove();
            mergeTasksInSession.add(traverser);
        }

//        logger.warn("Selected things to merge:");
//        printOutputMergeTasks(mergeTasksInSession);

        // Actually run through, merging the tasks in the working queue.
        for( ShardTraverser traverser: mergeTasksInSession ) {
            //logger.warn("*** Merging " + traverser.getIntervalsString());
            if( !traverser.isComplete() )
                traverser.waitForComplete();

            OutputMergeTask mergeTask = traverser.getOutputMergeTask();
            if( mergeTask != null ) {
                try {
                    mergeTask.merge();
                }
                catch(TribbleException ex) {
                    // Specifically catch Tribble I/O exceptions and rethrow them as Reviewed.  We don't expect
                    // any issues here because we created the Tribble output file mere moments ago and expect it to
                    // be completely valid.
                    throw new ReviewedGATKException("Unable to merge temporary Tribble output file.",ex);
                }
            }
        }

        long endTime = System.currentTimeMillis();

        totalOutputMergeTime += ( endTime - startTime );
    }

    /**
     * Queues the next traversal of a walker from the traversal tasks queue.
     *
     * @param walker     Walker to apply to the dataset.
     * @param reduceTree Tree of reduces to which to add this shard traverse.
     */
    protected void queueNextShardTraverse( Walker walker, ReduceTree reduceTree ) {
        if (!traversalTasks.hasNext())
            throw new IllegalStateException("Cannot traverse; no pending traversals exist.");

        final Shard shard = traversalTasks.next();

        // todo -- add ownership claim here

        final ShardTraverser traverser = new ShardTraverser(this, walker, shard, outputTracker);

        final Future traverseResult = threadPool.submit(traverser);

        // Add this traverse result to the reduce tree.  The reduce tree will call a callback to throw its entries on the queue.
        reduceTree.addEntry(traverseResult);
        outputMergeTasks.add(traverser);

//        logger.warn("adding merge task");
//        printOutputMergeTasks();

        // No more data?  Let the reduce tree know so it can finish processing what it's got.
        if (!isShardTraversePending())
            reduceTree.complete();
    }

    private synchronized void printOutputMergeTasks() {
        printOutputMergeTasks(outputMergeTasks);
    }

    private synchronized void printOutputMergeTasks(final Queue<ShardTraverser> tasks) {
        logger.info("Output merge tasks " + tasks.size());
        for ( final ShardTraverser traverser : tasks )
            logger.info(String.format("\t%s: complete? %b", traverser.getIntervalsString(), traverser.isComplete()));
    }

    /** Pulls the next reduce from the queue and runs it. */
    protected void queueNextTreeReduce( Walker walker ) {
        if (reduceTasks.size() == 0)
            throw new IllegalStateException("Cannot reduce; no pending reduces exist.");
        final TreeReduceTask reducer = reduceTasks.remove();
        reducer.setWalker((TreeReducible) walker);

        threadPool.submit(reducer);
    }

    /** Blocks until a free slot appears in the thread queue. */
    protected void waitForFreeQueueSlot() {
        final ThreadPoolMonitor monitor = new ThreadPoolMonitor();
        synchronized (monitor) {
            threadPool.submit(monitor);
            monitor.watch();
        }
    }

    /**
     * Callback for adding reduce tasks to the run queue.
     *
     * @return A new, composite future of the result of this reduce.
     */
    public Future notifyReduce( final Future lhs, final Future rhs ) {
        final TreeReduceTask reducer = new TreeReduceTask(new TreeReducer(this, lhs, rhs));
        reduceTasks.add(reducer);
        return reducer;
    }

    /**
     * Allows other threads to notify of an error during traversal.
     */
    protected synchronized RuntimeException notifyOfTraversalError(Throwable error) {
        return errorTracker.notifyOfError(error);
    }

    /** A small wrapper class that provides the TreeReducer interface along with the FutureTask semantics. */
    private class TreeReduceTask extends FutureTask {
        final private TreeReducer treeReducer;

        public TreeReduceTask( TreeReducer treeReducer ) {
            super(treeReducer);
            this.treeReducer = treeReducer;
        }

        public void setWalker( TreeReducible walker ) {
            treeReducer.setWalker(walker);
        }

        public boolean isReadyForReduce() {
            return treeReducer.isReadyForReduce();
        }
    }

    /**
     * Used by the ShardTraverser to report time consumed traversing a given shard.
     *
     * @param shardTraversalTime Elapsed time traversing a given shard.
     */
    synchronized void reportShardTraverseTime( long shardTraversalTime ) {
        totalShardTraverseTime += shardTraversalTime;
        totalCompletedTraversals++;
    }

    /**
     * Used by the TreeReducer to report time consumed reducing two shards.
     *
     * @param treeReduceTime Elapsed time reducing two shards.
     */
    synchronized void reportTreeReduceTime( long treeReduceTime ) {
        totalTreeReduceTime += treeReduceTime;
        totalCompletedTreeReduces++;

    }

    /** {@inheritDoc} */
    public int getNumberOfTasksInReduceQueue() {
        return reduceTasks.size();
    }

    /** {@inheritDoc} */
    public int getNumberOfTasksInIOQueue() {
        synchronized( outputMergeTasks ) {
            return outputMergeTasks.size();
        }
    }

    /** {@inheritDoc} */
    public long getTotalShardTraverseTimeMillis() {
        return totalShardTraverseTime;
    }

    /** {@inheritDoc} */
    public long getAvgShardTraverseTimeMillis() {
        if (totalCompletedTraversals == 0)
            return 0;
        return totalShardTraverseTime / totalCompletedTraversals;
    }

    /** {@inheritDoc} */
    public long getTotalTreeReduceTimeMillis() {
        return totalTreeReduceTime;
    }

    /** {@inheritDoc} */
    public long getAvgTreeReduceTimeMillis() {
        if (totalCompletedTreeReduces == 0)
            return 0;
        return totalTreeReduceTime / totalCompletedTreeReduces;
    }

    /** {@inheritDoc} */
    public long getTotalOutputMergeTimeMillis() {
        return totalOutputMergeTime;
    }
}
