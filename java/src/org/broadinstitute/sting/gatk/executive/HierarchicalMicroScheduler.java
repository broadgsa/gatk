package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.OutputTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.threading.ThreadPoolMonitor;

import java.io.File;
import java.util.List;
import java.util.Queue;
import java.util.LinkedList;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 5:41:04 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * A microscheduler that schedules shards according to a tree-like structure.
 * Requires a special walker tagged with a 'TreeReducible' interface.
 */
public class HierarchicalMicroScheduler extends MicroScheduler implements ReduceTree.TreeReduceNotifier {
    /**
     * How many outstanding output merges are allowed before the scheduler stops
     * allowing new processes and starts merging flat-out.
     */
    private static final int MAX_OUTSTANDING_OUTPUT_MERGES = 50;

    private TraverseLociByReference traversalEngine = null;

    /**
     * Manage currently running threads.
     */
    private ExecutorService threadPool;

    private Queue<Shard> traverseTasks = new LinkedList<Shard>();
    private Queue<TreeReduceTask> reduceTasks = new LinkedList<TreeReduceTask>();
    private Queue<OutputMerger> outputMergeTasks = new LinkedList<OutputMerger>();

    /**
     * Create a new hierarchical microscheduler to process the given reads and reference.
     * @param reads Reads file(s) to process.
     * @param refFile Reference for driving the traversal.
     * @param nThreadsToUse maximum number of threads to use to do the work
     */
    protected HierarchicalMicroScheduler( Walker walker, List<File> reads, File refFile, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods, int nThreadsToUse ) {
        super( walker, reads, refFile, rods );
        this.threadPool = Executors.newFixedThreadPool(nThreadsToUse);
        if( !(traversalEngine instanceof TraverseLociByReference) )
            throw new UnsupportedOperationException("Traversal engine supports only traverse loci by reference at this time.");
    }

    public void execute( Walker walker, List<GenomeLoc> intervals ) {
        // Fast fail for walkers not supporting TreeReducible interface.
        if( !(walker instanceof TreeReducible) )
            throw new IllegalArgumentException("Hierarchical microscheduler only works with TreeReducible walkers");

        ShardStrategy shardStrategy = getShardStrategy( reference, intervals );
        ReduceTree reduceTree = new ReduceTree( this );        

        walker.initialize();
        
        for(Shard shard: shardStrategy)
            traverseTasks.add(shard);

        while( isShardTraversePending() || isTreeReducePending() ) {
            // Too many files sitting around taking up space?  Merge them.
            if( isMergeLimitExceeded() )
                mergeExistingOutput();

            // Wait for the next slot in the queue to become free.
            waitForFreeQueueSlot();

            // Pick the next most appropriate task and run it.  In the interest of
            // memory conservation, hierarchical reduces always run before traversals.
            if( isTreeReduceReady() )
                queueNextTreeReduce( walker );
            else if( isShardTraversePending() )
                queueNextShardTraverse( walker, reduceTree );
        }

        // Merge any lingering output files.  If these files aren't ready,
        // sit around and wait for them, then merge them.
        mergeRemainingOutput();

        threadPool.shutdown();

        Object result = null;
        try {
            result = reduceTree.getResult().get();
        }
        catch(Exception ex) {
            throw new StingException("Unable to retrieve result", ex );
        }
        
        traversalEngine.printOnTraversalDone(result);
        walker.onTraversalDone(result);
    }

    /**
     * Returns true if there are unscheduled shard traversal waiting to run.
     * @return true if a shard traversal is waiting; false otherwise.
     */
    protected boolean isShardTraversePending() {
        return traverseTasks.size() > 0;
    }

    /**
     * Returns true if there are tree reduces that can be run without
     * blocking.
     * @return true if a tree reduce is ready; false otherwise.
     */
    protected boolean isTreeReduceReady() {
        if( reduceTasks.size() == 0 )
            return false;
        return reduceTasks.peek().isReadyForReduce();
    }

    /**
     * Returns true if there are tree reduces that need to be run before
     * the computation is complete.  Returns true if any entries are in the queue,
     * blocked or otherwise.
     * @return true if a tree reduce is pending; false otherwise.
     */
    protected boolean isTreeReducePending() {
        return reduceTasks.size() > 0;
    }

    /**
     * Returns whether the maximum number of files is sitting in the temp directory
     * waiting to be merged back in.
     * @return True if the merging needs to take priority.  False otherwise.
     */
    protected boolean isMergeLimitExceeded() {
        if( outputMergeTasks.size() < MAX_OUTSTANDING_OUTPUT_MERGES )
            return false;

        // If any of the first MAX_OUTSTANDING merges aren't ready, the merge limit
        // has not been exceeded.
        OutputMerger[] outputMergers = outputMergeTasks.toArray( new OutputMerger[0] );
        for( int i = 0; i < outputMergers.length; i++ ) {
            if( !outputMergers[i].isComplete() )
                return false;
        }

        // Everything's ready?  
        return true;
    }

    /**
     * Returns whether there is output waiting to be merged into the global output
     * streams right now.
     * @return True if this output is ready to be merged.  False otherwise.
     */
    protected boolean isOutputMergeReady() {
        return outputMergeTasks.peek().isComplete();
    }

    /**
     * Merging all output that's sitting ready in the OutputMerger queue into
     * the final data streams.
     */
    protected void mergeExistingOutput() {
        OutputTracker outputTracker = GenomeAnalysisEngine.instance.getOutputTracker();
        while( isOutputMergeReady() )
            outputMergeTasks.remove().mergeInto( outputTracker.getGlobalOutStream(), outputTracker.getGlobalErrStream() );
    }

    /**
     * Merge any output that hasn't yet been taken care of by the blocking thread.
     */
    protected void mergeRemainingOutput() {
        OutputTracker outputTracker = GenomeAnalysisEngine.instance.getOutputTracker();
        while( outputMergeTasks.size() > 0 ) {
            OutputMerger outputMerger = outputMergeTasks.remove();
            synchronized(outputMerger) {
                if( !outputMerger.isComplete() )
                    outputMerger.waitForOutputComplete();
            }
            outputMerger.mergeInto( outputTracker.getGlobalOutStream(), outputTracker.getGlobalErrStream() );            
        }
    }

    /**
     * Queues the next traversal of a walker from the traversal tasks queue.
     * @param walker Walker to apply to the dataset.
     * @param reduceTree Tree of reduces to which to add this shard traverse.
     */
    protected Future queueNextShardTraverse( Walker walker, ReduceTree reduceTree ) {
        if( traverseTasks.size() == 0 )
            throw new IllegalStateException( "Cannot traverse; no pending traversals exist.");

        Shard shard = traverseTasks.remove();
        OutputMerger outputMerger = new OutputMerger();

        ShardTraverser traverser = new ShardTraverser( traversalEngine,
                                                       walker,
                                                       shard,
                                                       getShardDataProvider(shard),
                                                       outputMerger );

        Future traverseResult = threadPool.submit(traverser);

        // Add this traverse result to the reduce tree.  The reduce tree will call a callback to throw its entries on the queue.
        reduceTree.addEntry( traverseResult );

        // No more data?  Let the reduce tree know so it can finish processing what it's got.
        if( !isShardTraversePending() )
            reduceTree.complete();

        outputMergeTasks.add(outputMerger);        

        return traverseResult;
    }

    /**
     * Pulls the next reduce from the queue and runs it.
     */
    protected void queueNextTreeReduce( Walker walker ) {
        if( reduceTasks.size() == 0 )
            throw new IllegalStateException( "Cannot reduce; no pending reduces exist.");
        TreeReduceTask reducer = reduceTasks.remove();
        reducer.setWalker( (TreeReducible)walker );

        threadPool.submit( reducer );
    }

    /**
     * Blocks until a free slot appears in the thread queue.
     */
    protected void waitForFreeQueueSlot() {
        ThreadPoolMonitor monitor = new ThreadPoolMonitor();
        synchronized(monitor) {
            threadPool.submit( monitor );
            monitor.watch();
        }
    }

    /**
     * Callback for adding reduce tasks to the run queue.
     * @return A new, composite future of the result of this reduce.
     */
    public Future notifyReduce( Future lhs, Future rhs ) {
        TreeReduceTask reducer = new TreeReduceTask( new TreeReducer( lhs, rhs ) );
        reduceTasks.add(reducer);
        return reducer; 
    }


    /**
     * A small wrapper class that provides the TreeReducer interface along with the FutureTask semantics.
     */
    private class TreeReduceTask extends FutureTask {
        private TreeReducer treeReducer = null;

        public TreeReduceTask( TreeReducer treeReducer ) {
            super(treeReducer);
            this.treeReducer = treeReducer;
        }

        public void setWalker( TreeReducible walker ) {
            treeReducer.setWalker( walker );
        }

        public boolean isReadyForReduce() {
            return treeReducer.isReadyForReduce();
        }
    }


}
