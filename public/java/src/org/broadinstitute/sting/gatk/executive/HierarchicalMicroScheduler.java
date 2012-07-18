package org.broadinstitute.sting.gatk.executive;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.ThreadLocalOutputTracker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.threading.ThreadPoolMonitor;

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
    private ThreadLocalOutputTracker outputTracker = new ThreadLocalOutputTracker();

    private final Queue<TreeReduceTask> reduceTasks = new LinkedList<TreeReduceTask>();

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
     * @param walker        the walker used to process the dataset.
     * @param reads         Reads file(s) to process.
     * @param reference     Reference for driving the traversal.
     * @param nThreadsToUse maximum number of threads to use to do the work
     */
    protected HierarchicalMicroScheduler(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, int nThreadsToUse ) {
        super(engine, walker, reads, reference, rods);
        this.threadPool = Executors.newFixedThreadPool(nThreadsToUse);
    }

    public Object execute( Walker walker, Iterable<Shard> shardStrategy ) {
        // Fast fail for walkers not supporting TreeReducible interface.
        if (!( walker instanceof TreeReducible ))
            throw new IllegalArgumentException("The GATK can currently run in parallel only with TreeReducible walkers");

        this.traversalTasks = shardStrategy.iterator();

        ReduceTree reduceTree = new ReduceTree(this);
        initializeWalker(walker);

        //
        // exception handling here is a bit complex.  We used to catch and rethrow exceptions all over
        // the place, but that just didn't work well.  Now we have a specific execution exception (inner class)
        // to use for multi-threading specific exceptions.  All RuntimeExceptions that occur within the threads are rethrown
        // up the stack as their underlying causes
        //
        while (isShardTraversePending() || isTreeReducePending()) {
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

        threadPool.shutdown();

        // Merge any lingering output files.  If these files aren't ready,
        // sit around and wait for them, then merge them.
        mergeExistingOutput(true);

        Object result = null;
        try {
            result = reduceTree.getResult().get();
            notifyTraversalDone(walker,result);
        }
        catch( InterruptedException ex ) { handleException(ex); }
        catch( ExecutionException ex ) { handleException(ex); }

        // do final cleanup operations
        outputTracker.close();
        cleanup();

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
            printOnTraversalDone(result);
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

        // Actually run through, merging the tasks in the working queue.
        for( ShardTraverser traverser: mergeTasksInSession ) {
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
                    final String reason = ex.getMessage();
                    throw new ReviewedStingException("Unable to merge temporary Tribble output file" + (reason == null ? "." : (" (" + reason + ").")), ex);
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

        Shard shard = traversalTasks.next();

        // todo -- add ownership claim here

        ShardTraverser traverser = new ShardTraverser(this,
                traversalEngine,
                walker,
                shard,
                outputTracker);

        Future traverseResult = threadPool.submit(traverser);

        // Add this traverse result to the reduce tree.  The reduce tree will call a callback to throw its entries on the queue.
        reduceTree.addEntry(traverseResult);
        outputMergeTasks.add(traverser);

        // No more data?  Let the reduce tree know so it can finish processing what it's got.
        if (!isShardTraversePending())
            reduceTree.complete();
    }

    /** Pulls the next reduce from the queue and runs it. */
    protected void queueNextTreeReduce( Walker walker ) {
        if (reduceTasks.size() == 0)
            throw new IllegalStateException("Cannot reduce; no pending reduces exist.");
        TreeReduceTask reducer = reduceTasks.remove();
        reducer.setWalker((TreeReducible) walker);

        threadPool.submit(reducer);
    }

    /** Blocks until a free slot appears in the thread queue. */
    protected void waitForFreeQueueSlot() {
        ThreadPoolMonitor monitor = new ThreadPoolMonitor();
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
    public Future notifyReduce( Future lhs, Future rhs ) {
        TreeReduceTask reducer = new TreeReduceTask(new TreeReducer(this, lhs, rhs));
        reduceTasks.add(reducer);
        return reducer;
    }

    /**
     * Handle an exception that occurred in a worker thread as needed by this scheduler.
     *
     * The way to use this function in a worker is:
     *
     * try { doSomeWork();
     * catch ( InterruptedException ex ) { hms.handleException(ex); }
     * catch ( ExecutionException ex ) { hms.handleException(ex); }
     *
     * @param ex the exception that occurred in the worker thread
     */
    protected final void handleException(InterruptedException ex) {
        throw new HierarchicalMicroScheduler.ExecutionFailure("Hierarchical reduce interrupted", ex);
    }

    /**
     * Handle an exception that occurred in a worker thread as needed by this scheduler.
     *
     * The way to use this function in a worker is:
     *
     * try { doSomeWork();
     * catch ( InterruptedException ex ) { hms.handleException(ex); }
     * catch ( ExecutionException ex ) { hms.handleException(ex); }
     *
     * @param ex the exception that occurred in the worker thread
     */
    protected final void handleException(ExecutionException ex) {
        if ( ex.getCause() instanceof RuntimeException )
            // if the cause was a runtime exception that's what we want to send up the stack
            throw (RuntimeException )ex.getCause();
        else
            throw new HierarchicalMicroScheduler.ExecutionFailure("Hierarchical reduce failed", ex);
    }



    /** A small wrapper class that provides the TreeReducer interface along with the FutureTask semantics. */
    private class TreeReduceTask extends FutureTask {
        private TreeReducer treeReducer = null;

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
     * A specific exception class for HMS-specific failures such as
     * Interrupted or ExecutionFailures that aren't clearly the fault
     * of the underlying walker code
     */
    public static class ExecutionFailure extends ReviewedStingException {
        public ExecutionFailure(final String s, final Throwable throwable) {
            super(s, throwable);
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
