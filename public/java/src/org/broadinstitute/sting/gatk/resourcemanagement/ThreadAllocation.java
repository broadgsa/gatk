/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.resourcemanagement;

import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Models how threads are distributed between various components of the GATK.
 */
public class ThreadAllocation {
    /**
     * The number of CPU threads to be used by the GATK.
     */
    private final int numCPUThreads;

    /**
     * Number of threads to devote exclusively to IO.  Default is 0.
     */
    private final int numIOThreads;

    /**
     * Should we monitor thread efficiency?
     */
    private final boolean monitorThreads;

    public int getNumCPUThreads() {
        return numCPUThreads;
    }

    public int getNumIOThreads() {
        return numIOThreads;
    }

    public boolean shouldMonitorThreads() {
        return monitorThreads;
    }

    /**
     * Construct the default thread allocation.
     */
    public ThreadAllocation() {
        this(1, null, null, false);
    }

    /**
     * Set up the thread allocation.  Default allocation is 1 CPU thread, 0 IO threads.
     * (0 IO threads means that no threads are devoted exclusively to IO; they're inline on the CPU thread).
     * @param totalThreads Complete number of threads to allocate.
     * @param numCPUThreads Total number of threads allocated to the traversal.
     * @param numIOThreads Total number of threads allocated exclusively to IO.
     */
    public ThreadAllocation(final int totalThreads, final Integer numCPUThreads, final Integer numIOThreads, final boolean monitorThreads) {
        // If no allocation information is present, allocate all threads to CPU
        if(numCPUThreads == null && numIOThreads == null) {
            this.numCPUThreads = totalThreads;
            this.numIOThreads = 0;
        }
        // If only CPU threads are specified, allocate remainder to IO (minimum 0 dedicated IO threads).
        else if(numIOThreads == null) {
            if(numCPUThreads > totalThreads)
                throw new UserException(String.format("Invalid thread allocation.  User requested %d threads in total, but the count of cpu threads (%d) is higher than the total threads",totalThreads,numCPUThreads));
            this.numCPUThreads = numCPUThreads;
            this.numIOThreads = totalThreads - numCPUThreads;
        }
        // If only IO threads are specified, allocate remainder to CPU (minimum 1 dedicated CPU thread).
        else if(numCPUThreads == null) {
            if(numIOThreads > totalThreads)
                throw new UserException(String.format("Invalid thread allocation.  User requested %d threads in total, but the count of io threads (%d) is higher than the total threads",totalThreads,numIOThreads));
            this.numCPUThreads = Math.max(1,totalThreads-numIOThreads);
            this.numIOThreads = numIOThreads;
        }
        else {
            if(numCPUThreads + numIOThreads != totalThreads)
                throw new UserException(String.format("Invalid thread allocation.  User requested %d threads in total, but the count of cpu threads (%d) + the count of io threads (%d) does not match",totalThreads,numCPUThreads,numIOThreads));
            this.numCPUThreads = numCPUThreads;
            this.numIOThreads = numIOThreads;
        }

        this.monitorThreads = monitorThreads;
    }
}
