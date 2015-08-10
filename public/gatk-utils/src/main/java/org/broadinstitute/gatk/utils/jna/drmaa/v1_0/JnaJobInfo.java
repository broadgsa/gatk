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

package org.broadinstitute.gatk.utils.jna.drmaa.v1_0;

import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.JobInfo;

import java.util.Map;

/**
 * JNA mapping from Java to C DRMAA binding.
 */
public class JnaJobInfo implements JobInfo {

    private final String jobId;
    private final Map<String, String> rusage;
    private final boolean hasExited;
    private final int exitStatus;
    private final boolean hasSignaled;
    private final String terminatingSignal;
    private final boolean hasCoreDump;
    private final boolean wasAborted;
            
    public JnaJobInfo(String jobId, Map<String, String> rusage, boolean hasExited, int exitStatus, boolean hasSignaled, String terminatingSignal, boolean hasCoreDump, boolean wasAborted) {
        this.jobId = jobId;
        this.rusage = rusage;
        this.hasExited = hasExited;
        this.exitStatus = exitStatus;
        this.hasSignaled = hasSignaled;
        this.terminatingSignal = terminatingSignal;
        this.hasCoreDump = hasCoreDump;
        this.wasAborted = wasAborted;
    }

    @Override
    public String getJobId() throws DrmaaException {
        return this.jobId;
    }

    @Override
    public Map getResourceUsage() throws DrmaaException {
        return rusage;
    }

    @Override
    public boolean hasExited() throws DrmaaException {
        return hasExited;
    }

    @Override
    public int getExitStatus() throws DrmaaException {
        if (!hasExited)
            throw new IllegalStateException("job has not exited");
        return exitStatus;
    }

    @Override
    public boolean hasSignaled() throws DrmaaException {
        return hasSignaled;
    }

    @Override
    public String getTerminatingSignal() throws DrmaaException {
        if (!hasSignaled)
            throw new IllegalStateException("job has not signaled");
        return terminatingSignal;
    }

    @Override
    public boolean hasCoreDump() throws DrmaaException {
        return hasCoreDump;
    }

    @Override
    public boolean wasAborted() throws DrmaaException {
        return wasAborted;
    }
}
