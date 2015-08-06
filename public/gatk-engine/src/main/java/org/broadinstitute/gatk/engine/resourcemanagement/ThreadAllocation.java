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

package org.broadinstitute.gatk.engine.resourcemanagement;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

/**
 * Models how threads are distributed between various components of the GATK.
 */
public class ThreadAllocation {
    /**
     * The number of CPU threads to be used by the GATK.
     */
    private final int numDataThreads;

    /**
     * The number of CPU threads per data thread for GATK processing
     */
    private final int numCPUThreadsPerDataThread;

    /**
     * Number of threads to devote exclusively to IO.  Default is 0.
     */
    private final int numIOThreads;

    /**
     * Should we monitor thread efficiency?
     */
    private final boolean monitorEfficiency;

    public int getNumDataThreads() {
        return numDataThreads;
    }

    public int getNumCPUThreadsPerDataThread() {
        return numCPUThreadsPerDataThread;
    }

    public int getNumIOThreads() {
        return numIOThreads;
    }

    public boolean monitorThreadEfficiency() {
        return monitorEfficiency;
    }

    /**
     * Are we running in parallel mode?
     *
     * @return true if any parallel processing is enabled
     */
    public boolean isRunningInParallelMode() {
        return getTotalNumThreads() > 1;
    }

    /**
     * What is the total number of threads in use by the GATK?
     *
     * @return the sum of all thread allocations in this object
     */
    public int getTotalNumThreads() {
        return getNumDataThreads() * getNumCPUThreadsPerDataThread() + getNumIOThreads();
    }

    /**
     * Construct the default thread allocation.
     */
    public ThreadAllocation() {
        this(1, 1, 0, false);
    }

    /**
     * Set up the thread allocation.  Default allocation is 1 CPU thread, 0 IO threads.
     * (0 IO threads means that no threads are devoted exclusively to IO; they're inline on the CPU thread).
     * @param numDataThreads Total number of threads allocated to the traversal.
     * @param numCPUThreadsPerDataThread The number of CPU threads per data thread to allocate
     * @param numIOThreads Total number of threads allocated exclusively to IO.
     * @param monitorEfficiency should we monitor threading efficiency in the GATK?
     */
    public ThreadAllocation(final int numDataThreads,
                            final int numCPUThreadsPerDataThread,
                            final int numIOThreads,
                            final boolean monitorEfficiency) {
        if ( numDataThreads < 1 ) throw new ReviewedGATKException("numDataThreads cannot be less than 1, but saw " + numDataThreads);
        if ( numCPUThreadsPerDataThread < 1 ) throw new ReviewedGATKException("numCPUThreadsPerDataThread cannot be less than 1, but saw " + numCPUThreadsPerDataThread);
        if ( numIOThreads < 0 ) throw new ReviewedGATKException("numIOThreads cannot be less than 0, but saw " + numIOThreads);

        this.numDataThreads = numDataThreads;
        this.numCPUThreadsPerDataThread = numCPUThreadsPerDataThread;
        this.numIOThreads = numIOThreads;
        this.monitorEfficiency = monitorEfficiency;
    }
}
