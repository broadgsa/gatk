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

package org.broadinstitute.gatk.utils.nanoScheduler;

/**
 * Holds the results of a map job suitable for producer/consumer threading
 * via a BlockingQueue
 */
class MapResult<MapType> extends EOFMarkedValue<MapType> implements Comparable<MapResult<MapType>> {
    final int jobID;

    /**
     * Create a new MapResult with value datum and jod jobID ID
     *
     * @param datum the value produced by the map job
     * @param jobID the id of the map job (for correctness testing)
     */
    MapResult(final MapType datum, final int jobID) {
        super(datum);
        this.jobID = jobID;
        if ( jobID < 0 ) throw new IllegalArgumentException("JobID must be >= 0");
    }

    MapResult(final int jobID) {
        super();
        this.jobID = jobID;
        if ( jobID < 0 ) throw new IllegalArgumentException("JobID must be >= 0");
    }

    /**
     * @return the job ID of the map job that produced this MapResult
     */
    public int getJobID() {
        return jobID;
    }

    /**
     * Compare these MapResults in order of JobID.
     *
     * @param o
     * @return
     */
    @Override
    public int compareTo(MapResult<MapType> o) {
        return Integer.valueOf(jobID).compareTo(o.getJobID());
    }

    @Override
    public String toString() {
        return "[MapResult id=" + jobID + "]";
    }
}
