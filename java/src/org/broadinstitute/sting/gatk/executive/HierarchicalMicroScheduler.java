package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.List;
import java.util.Queue;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.ExecutionException;
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

        while( isShardTraversePending() || isTreeReducePending() ) {
            waitForFreeQueueSlot();

            if( isTreeReduceReady() )
                queueNextTreeReduce( walker );
            else {
                Future traverseResult = queueNextShardTraverse( walker, dataSource );

                // Add this traverse result to the reduce tree.  The reduce tree will call a callback to throw its entries on the queue.
                reduceTree.addEntry( traverseResult );

                // No more data?  Let the reduce tree know so it can finish processing what it's got.
                if( !isShardTraversePending() )
                    reduceTree.complete();
            }
        }

        Object result = reduceTree.getResult();
        
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
     * Queues the next traversal of a walker from the traversal tasks queue.
     * @param walker Walker to apply to the dataset.
     * @param dataSource Source of the reads
     */
    protected Future queueNextShardTraverse( Walker walker, SAMDataSource dataSource ) {
        if( traverseTasks.size() == 0 )
            throw new IllegalStateException( "Cannot traverse; no pending traversals exist.");
        ShardTraverser traverser = new ShardTraverser( walker, traverseTasks.remove(), dataSource );
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
     * Carries the walker over a given shard.
     */
    private class ShardTraverser implements Callable {
        private Walker walker;
        private Shard shard;
        private SAMDataSource dataSource;

        public ShardTraverser( Walker walker, Shard shard, SAMDataSource dataSource ) {
            this.walker = walker;
            this.shard = shard;
            this.dataSource = dataSource;            
        }

        public Object call() {
            GenomeLoc span = shard.getGenomeLoc();
            Object accumulator = ((LocusWalker<?,?>)walker).reduceInit();

            MergingSamRecordIterator2 readShard = null;
            try {
                readShard = dataSource.seek( span );
            }
            catch( SimpleDataSourceLoadException ex ) {
                throw new RuntimeException( ex );
            }

            ReferenceProvider referenceProvider = new ReferenceProvider( reference, span );
            LocusContextProvider locusProvider = new LocusContextProvider( readShard );

            accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );

            readShard.close();

            return accumulator;
        }
    }

    /**
     * Waits for a signal to come through that the thread pool has run
     * a given task and therefore has a free slot.
     */
    private class ThreadPoolMonitor implements Runnable {
        public synchronized void watch() {
            try {
                wait();
            }
            catch( InterruptedException ex ) {
                logger.error("ThreadPoolMonitor interrupted:" + ex.getStackTrace());
                throw new RuntimeException("ThreadPoolMonitor interrupted", ex);
            }
        }

        public synchronized void run() {
            notify();
        }
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


    /**
     * Represents a 'potential' reduce...a reduce that will be ready at some point in the future.
     * Provides services for indicating when all data is prepared for the reduce and services to
     * actually make that reduce happen.
     */
    private class TreeReducer implements Callable {
        private TreeReducible walker;
        private final Future lhs;
        private final Future rhs;

        public TreeReducer( Future lhs ) {
            this( lhs, null );
        }

        public TreeReducer( Future lhs, Future rhs ) {
            this.lhs = lhs;
            this.rhs = rhs;
        }

        public void setWalker( TreeReducible walker ) {
            this.walker = walker;
        }

        public boolean isReadyForReduce() {
            if( lhs == null )
                throw new IllegalStateException(String.format("Insufficient data on which to reduce; lhs = %s, rhs = %s", lhs, rhs) );

            if( rhs == null )
                return lhs.isDone();

            return lhs.isDone() && rhs.isDone();
        }

        public Object call() {
            try {
                if( lhs == null )
                    return lhs.get();
                else
                    return walker.reduce( lhs.get(), rhs.get() );
            }
            catch( InterruptedException ex ) {
                throw new RuntimeException("Hierarchical reduce interrupted", ex);
            }
            catch( ExecutionException ex ) {
                throw new RuntimeException("Hierarchical reduce failed", ex);
            }
        }
    }

}
