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

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.providers.ShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.io.ThreadGroupOutputTracker;
import org.broadinstitute.gatk.engine.traversals.TraversalEngine;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
    final private HierarchicalMicroScheduler microScheduler;
    final private Walker walker;
    final private Shard shard;
    final private ThreadGroupOutputTracker outputTracker;
    private OutputMergeTask outputMergeTask;

    /** our log, which we want to capture anything from this class */
    final protected static Logger logger = Logger.getLogger(ShardTraverser.class);

    /**
     * Is this traversal complete?
     */
    private boolean complete = false;

    public ShardTraverser( HierarchicalMicroScheduler microScheduler,
                           Walker walker,
                           Shard shard,
                           ThreadGroupOutputTracker outputTracker) {
        this.microScheduler = microScheduler;
        this.walker = walker;
        this.shard = shard;
        this.outputTracker = outputTracker;
    }

    public Object call() {
        final Object traversalEngineKey = Thread.currentThread();
        final TraversalEngine traversalEngine = microScheduler.borrowTraversalEngine(traversalEngineKey);

        try {
            final long startTime = System.currentTimeMillis();

            // this is CRITICAL -- initializes output maps in this master thread,
            // so that any subthreads created by the traversal itself can access this map
            outputTracker.initializeStorage();

            Object accumulator = walker.reduceInit();
            final WindowMaker windowMaker = new WindowMaker(shard,microScheduler.getEngine().getGenomeLocParser(),
                    microScheduler.getReadIterator(shard),
                    shard.getGenomeLocs(),
                    microScheduler.engine.getSampleDB().getSampleNames()); // todo: microScheduler.engine is protected - is it okay to user it here?

            for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
                final ShardDataProvider dataProvider = new LocusShardDataProvider(shard,iterator.getSourceInfo(),microScheduler.getEngine().getGenomeLocParser(),iterator.getLocus(),iterator,microScheduler.reference,microScheduler.rods);
                accumulator = traversalEngine.traverse(walker, dataProvider, accumulator);
                dataProvider.close();
            }

            windowMaker.close();
            outputMergeTask = outputTracker.closeStorage();

            final long endTime = System.currentTimeMillis();

            microScheduler.reportShardTraverseTime(endTime-startTime);

            return accumulator;
        } catch(Throwable t) {
            // Notify that an exception has occurred and rethrow it.
            throw microScheduler.notifyOfTraversalError(t);
        } finally {
            synchronized(this) {
                complete = true;
                microScheduler.returnTraversalEngine(traversalEngineKey, traversalEngine);
                notifyAll();
            }
        }
    }

    /**
     * Return a human readable string describing the intervals this traverser is operating on
     * @return
     */
    public String getIntervalsString() {
        return Utils.join(",", shard.getGenomeLocs());
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
            throw new ReviewedGATKException("Interrupted while waiting for more output to be finalized.",ex);
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
