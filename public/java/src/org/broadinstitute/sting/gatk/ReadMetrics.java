/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk;

import net.sf.picard.filter.SamRecordFilter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Collections;
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
    private long nSkippedReads;
    private long nUnmappedReads;
    private long nNotPrimary;
    private long nBadAlignments;
    private long nSkippedIndels;
    private long nDuplicates;
    private Map<Class, Long> counter = new HashMap<Class, Long>();

    /**
     * Combines these metrics with a set of other metrics, storing the results in this class.
     * @param metrics The metrics to fold into this class.
     */
    public void incrementMetrics(ReadMetrics metrics) {
        nRecords += metrics.nRecords;
        nReads += metrics.nReads;
        nSkippedReads += metrics.nSkippedReads;
        nUnmappedReads += metrics.nUnmappedReads;
        nNotPrimary += metrics.nNotPrimary;
        nBadAlignments += metrics.nBadAlignments;
        nSkippedIndels += metrics.nSkippedIndels;
        nDuplicates += metrics.nDuplicates;
        for(Map.Entry<Class,Long> counterEntry: metrics.counter.entrySet()) {
            Class counterType = counterEntry.getKey();
            long newValue = (counter.containsKey(counterType) ? counter.get(counterType) : 0) + counterEntry.getValue();
            counter.put(counterType,newValue);
        }
    }

    /**
     * Create a copy of the given read metrics.
     * @return
     */
    public ReadMetrics clone() {
        ReadMetrics newMetrics;
        try {
            newMetrics = (ReadMetrics)super.clone();
        }
        catch(CloneNotSupportedException ex) {
            throw new ReviewedStingException("Unable to clone runtime metrics",ex);
        }
        newMetrics.nRecords = nRecords;
        newMetrics.nReads = nReads;
        newMetrics.nSkippedReads = nSkippedReads;
        newMetrics.nUnmappedReads = nUnmappedReads;
        newMetrics.nNotPrimary = nNotPrimary;
        newMetrics.nBadAlignments = nBadAlignments;
        newMetrics.nSkippedIndels = nSkippedIndels;
        newMetrics.nDuplicates = nDuplicates;
        newMetrics.counter = new HashMap<Class,Long>(counter);

        return newMetrics;
    }


    public void incrementFilter(SamRecordFilter filter) {
        long c = 0;
        if ( counter.containsKey(filter.getClass()) ) {
            c = counter.get(filter.getClass());
        }

        counter.put(filter.getClass(), c + 1L);
    }

    public Map<String,Long> getCountsByFilter() {
        final TreeMap<String, Long> sortedCounts = new TreeMap<String, Long>();
        for(Map.Entry<Class,Long> counterEntry: counter.entrySet()) {
            sortedCounts.put(counterEntry.getKey().getSimpleName(),counterEntry.getValue());
        }
        return sortedCounts;
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
    public void incrementNumIterations() {
        nRecords++;
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

    /**
     * Gets the cumulative number of reads skipped in the course of this run.
     * @return Cumulative number of reads skipped in the course of this run.
     */
    public long getNumSkippedReads() {
        return nSkippedReads;
    }

    /**
     * Increments the cumulative number of reads skipped in the course of this run.
     */
    public void incrementNumSkippedReads() {
        nSkippedReads++;
    }

    /**
     * Gets the number of unmapped reads skipped in the course of this run.
     * @return The number of unmapped reads skipped.
     */
    public long getNumUnmappedReads() {
        return nUnmappedReads;
    }

    /**
     * Increments the number of unmapped reads skipped in the course of this run.
     */
    public void incrementNumUnmappedReads() {
        nUnmappedReads++;
    }

    /**
     *
     * @return
     */
    public long getNumNonPrimaryReads() {
        return nNotPrimary;
    }

    /**
     *
     */
    public void incrementNumNonPrimaryReads() {
        nNotPrimary++;
    }

    /**
     *
     * @return
     */
    public long getNumBadAlignments() {
        return nBadAlignments;
    }

    /**
     *
     */
    public void incrementNumBadAlignments() {
        nBadAlignments++;
    }

    /**
     *
     * @return
     */
    public long getNumSkippedIndels() {
        return nSkippedIndels;
    }

    /**
     *
     */
    public void incrementNumSkippedIndels() {
        nSkippedIndels++;
    }

    /**
     *
     * @return
     */
    public long getNumDuplicates() {
        return nDuplicates;
    }

    /**
     *
     */
    public void incrementNumDuplicates() {
        nDuplicates++;
    }

}
