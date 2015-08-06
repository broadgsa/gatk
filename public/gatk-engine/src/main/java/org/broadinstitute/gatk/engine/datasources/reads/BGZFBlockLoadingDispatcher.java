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

package org.broadinstitute.gatk.engine.datasources.reads;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.LinkedList;
import java.util.Queue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Preloads BGZF blocks in preparation for unzipping and data processing.
 * TODO: Right now, the block loader has all threads blocked waiting for a work request.  Ultimately this should
 * TODO: be replaced with a central thread management strategy.
 */
public class BGZFBlockLoadingDispatcher {
    /**
     * The file handle cache, used when allocating blocks from the dispatcher.
     */
    private final FileHandleCache fileHandleCache;

    private final ExecutorService threadPool;

    private final Queue<BAMAccessPlan> inputQueue;

    public BGZFBlockLoadingDispatcher(final int numThreads, final int numFileHandles) {
        threadPool = Executors.newFixedThreadPool(numThreads);
        fileHandleCache = new FileHandleCache(numFileHandles);
        inputQueue = new LinkedList<BAMAccessPlan>();

        threadPool.execute(new BlockLoader(this,fileHandleCache,true));
    }

    /**
     * Initiates a request for a new block load.
      * @param readerPosition Position at which to load.
     */
    void queueBlockLoad(final BAMAccessPlan readerPosition) {
        synchronized(inputQueue) {
            inputQueue.add(readerPosition);
            inputQueue.notify();
        }
    }

    /**
     * Claims the next work request from the queue.
     * @return The next work request, or null if none is available.
     */
    BAMAccessPlan claimNextWorkRequest() {
        synchronized(inputQueue) {
            while(inputQueue.isEmpty()) {
                try {
                    inputQueue.wait();
                }
                catch(InterruptedException ex) {
                    throw new ReviewedGATKException("Interrupt occurred waiting for next block reader work item");
                }
            }
            return inputQueue.poll();
        }
    }
}
