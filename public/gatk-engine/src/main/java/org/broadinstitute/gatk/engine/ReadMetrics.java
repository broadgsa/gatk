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

package org.broadinstitute.gatk.engine;

import htsjdk.samtools.filter.SamRecordFilter;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Holds a bunch of basic information about the traversal.
 */
public class ReadMetrics implements Cloneable {
    // Number of records (loci, reads) we've processed
    private long nRecords;
    // How many reads have we processed, along with those skipped for various reasons
    private long nReads;

    // keep track of filtered records by filter type (class)
    private Map<String, Long> filterCounter = new HashMap<>();

    /**
     * Combines these metrics with a set of other metrics, storing the results in this class.
     * @param metrics The metrics to fold into this class.
     */
    public synchronized void incrementMetrics(ReadMetrics metrics) {
        nRecords += metrics.nRecords;
        nReads += metrics.nReads;
        for(Map.Entry<String, Long> counterEntry: metrics.filterCounter.entrySet()) {
            final String counterType = counterEntry.getKey();
            final long newValue = (filterCounter.containsKey(counterType) ? filterCounter.get(counterType) : 0) + counterEntry.getValue();
            filterCounter.put(counterType, newValue);
        }
    }

    /**
     * Create a copy of the given read metrics.
     * @return a non-null clone
     */
    public ReadMetrics clone() {
        ReadMetrics newMetrics;
        try {
            newMetrics = (ReadMetrics)super.clone();
        }
        catch(CloneNotSupportedException ex) {
            throw new ReviewedGATKException("Unable to clone runtime metrics",ex);
        }
        newMetrics.nRecords = nRecords;
        newMetrics.nReads = nReads;
        newMetrics.filterCounter = new HashMap<>(filterCounter);

        return newMetrics;
    }


    public void setFilterCount(final String filter, final long count) {
        filterCounter.put(filter, count);
    }

    public Map<String,Long> getCountsByFilter() {
        return new TreeMap<>(filterCounter);
    }

    /**
     * Gets the number of 'iterations' (one call of filter/map/reduce sequence) performed.
     * @return The number of iterations completed.
     */
    public long getNumIterations() {
        return nRecords;
    }

    /**
     * Increments the number of 'iterations' (one call of filter/map/reduce sequence) completed.
     */
    public void incrementNumIterations(final long by) {
        nRecords += by;
    }

    /**
     * Increments the number of 'iterations' (one call of filter/map/reduce sequence) completed.
     */
    public void incrementNumIterations() {
        incrementNumIterations(1);
    }

    public long getNumReadsSeen() {
        return nReads;
    }

    /**
     * Increments the number of reads seen in the course of this run.
     */
    public void incrementNumReadsSeen() {
        nReads++;
    }
}
