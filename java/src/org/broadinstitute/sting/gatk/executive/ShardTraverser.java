package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.io.ThreadLocalOutputTracker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.StingException;

import java.util.concurrent.Callable;
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 4:40:38 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */
/**
 * Carries the walker over a given shard, in a callable interface.
 */
public class ShardTraverser implements Callable {
    private HierarchicalMicroScheduler microScheduler;
    private Walker walker;
    private TraversalEngine traversalEngine;
    private ShardDataProvider dataProvider;
    private ThreadLocalOutputTracker outputTracker;
    private OutputMergeTask outputMergeTask;

    /**
     * Is this traversal complete?
     */
    private boolean complete = false;

    public ShardTraverser( HierarchicalMicroScheduler microScheduler,
                           TraversalEngine traversalEngine,
                           Walker walker,
                           ShardDataProvider dataProvider,
                           ThreadLocalOutputTracker outputTracker ) {
        this.microScheduler = microScheduler;
        this.walker = walker;
        this.traversalEngine = traversalEngine;
        this.dataProvider = dataProvider;
        this.outputTracker = outputTracker;
    }

    public Object call() {
        long startTime = System.currentTimeMillis(); 

        Object accumulator = walker.reduceInit();
        try {
            accumulator = traversalEngine.traverse( walker, dataProvider, accumulator );
        }
        finally {
            dataProvider.close();
            outputMergeTask = outputTracker.closeStorage();

            synchronized(this) {
                complete = true;
                notifyAll();
            }
        }

        long endTime = System.currentTimeMillis();

        microScheduler.reportShardTraverseTime(endTime-startTime);

        return accumulator;
    }

    /**
     * Has this traversal completed?
     * @return True if completed, false otherwise.
     */
    public boolean isComplete() {
        synchronized(this) {
            return complete;
        }
    }

   /**
     * Waits for any the given OutputMerger to be ready for merging.
     */
    public void waitForComplete() {
        try {
            synchronized(this) {
                if( isComplete() )
                    return;
                wait();
            }
        }
        catch( InterruptedException ex ) {
            throw new StingException("Interrupted while waiting for more output to be finalized.",ex);
        }
    }

    /**
     * Gets the output merge task associated with the given shard.
     * @return OutputMergeTask if one exists; null if nothing needs to be merged.
     */
    public OutputMergeTask getOutputMergeTask() {
        return outputMergeTask;
    }
}
