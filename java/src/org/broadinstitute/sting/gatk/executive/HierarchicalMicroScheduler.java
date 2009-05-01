package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.GenomeAnalysisTK;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.threading.ThreadPoolMonitor;

import java.io.File;
import java.io.OutputStream;
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
    protected HierarchicalMicroScheduler( List<File> reads, File refFile, int nThreadsToUse ) {
        super( reads, refFile );


        this.threadPool = Executors.newFixedThreadPool(nThreadsToUse);        
        traversalEngine = new TraverseLociByReference( reads, refFile, new java.util.ArrayList() );
    }

    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    public void execute( Walker walker, List<GenomeLoc> intervals ) {
        // Fast fail for walkers not supporting TreeReducible interface.
        if( !(walker instanceof TreeReducible) )
            throw new IllegalArgumentException("Hierarchical microscheduler only works with TreeReducible walkers");

        ShardStrategy shardStrategy = getShardStrategy( reference, intervals );
        SAMDataSource dataSource = getReadsDataSource();

        ReduceTree reduceTree = new ReduceTree( this );        

        walker.initialize();
        
        for(Shard shard: shardStrategy)
            traverseTasks.add(shard);

        while( isShardTraversePending() || isTreeReducePending() || isOutputMergePending() ) {
            waitForFreeQueueSlot();

            if( isTreeReduceReady() )
                queueNextTreeReduce( walker );
            else if( isShardTraversePending() ) {
                Future traverseResult = queueNextShardTraverse( walker, dataSource );

                // Add this traverse result to the reduce tree.  The reduce tree will call a callback to throw its entries on the queue.
                reduceTree.addEntry( traverseResult );

                // No more data?  Let the reduce tree know so it can finish processing what it's got.
                if( !isShardTraversePending() )
                    reduceTree.complete();
            }
            else if( isOutputMergeReady() ) {
                queueNextOutputMerge();
            }
        }

        threadPool.shutdown();

        Object result = null;
        try {
            result = reduceTree.getResult().get();
        }
        catch(Exception ex) {
            throw new StingException("Unable to retrieve result", ex );
        }
        
        traversalEngine.printOnTraversalDone("loci", result);
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
     * Returns whether there is output waiting to be merged into the global output
     * streams right now.
     * @return True if this output is ready to be merged.  False otherwise.
     */
    protected boolean isOutputMergeReady() {
        return !OutputMerger.isMergeQueued() && outputMergeTasks.peek().isComplete();
    }

    /**
     * Returns whether there is output that will eventually need to be merged into
     * the output streams.
     * @return True if output will eventually need to be merged.  False otherwise.
     */
    protected boolean isOutputMergePending() {
        return outputMergeTasks.size() > 0;
    }

    /**
     * Queues the next traversal of a walker from the traversal tasks queue.
     * @param walker Walker to apply to the dataset.
     * @param dataSource Source of the reads
     */
    protected Future queueNextShardTraverse( Walker walker, SAMDataSource dataSource ) {
        if( traverseTasks.size() == 0 )
            throw new IllegalStateException( "Cannot traverse; no pending traversals exist.");

        ShardOutput shardOutput = new ShardOutput();

        ShardTraverser traverser = new ShardTraverser( traversalEngine,
                                                       walker,
                                                       traverseTasks.remove(),
                                                       reference,
                                                       dataSource,
                                                       shardOutput );

        outputMergeTasks.add(new OutputMerger(shardOutput,
                                              GenomeAnalysisTK.Instance.getOutputTracker().getGlobalOutStream(),
                                              GenomeAnalysisTK.Instance.getOutputTracker().getGlobalErrStream()));

        return threadPool.submit(traverser);
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
     * Pulls the next output merge and puts it on the queue.
     */
    protected void queueNextOutputMerge() {
        if( outputMergeTasks.size() == 0 )
            throw new IllegalStateException( "Cannot merge output; no pending merges exist.");
        if( OutputMerger.isMergeQueued() )
            throw new IllegalStateException( "Cannot merge output; another merge has already been queued.");

        OutputMerger.queueMerge();
        threadPool.submit( outputMergeTasks.remove() );
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
