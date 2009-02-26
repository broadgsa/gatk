/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.filter;

import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.picard.util.CloserUtil;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Filtering Iterator which takes a filter and an iterator and iterates
 * through only those records which are not rejected by the filter.
 *
 * @author Kathleen Tibbetts
 */
public class FilteringIterator implements CloseableIterator<SAMRecord> {

    private final Iterator<SAMRecord> iterator;
    private final SamRecordFilter filter;
    private SAMRecord next = null;

    /**
     * Constructor
     *
     * @param iterator  the backing iterator
     * @param filter    the filter (which may be a FilterAggregator)
     */
    public FilteringIterator(Iterator<SAMRecord> iterator, SamRecordFilter filter) {
        this.iterator = iterator;
        this.filter = filter;
        next = getNextRecord();
    }

    /**
     * Returns true if the iteration has more elements.
     *
     * @return  true if the iteration has more elements.  Otherwise returns false.
     */
    public boolean hasNext() {
        return next != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return  the next element in the iteration
     * @throws java.util.NoSuchElementException
     */
    public SAMRecord next() {
        if (next == null) {
            throw new NoSuchElementException("Iterator has no more elements.");
        }
        SAMRecord result = next;
        next = getNextRecord();
        return result;
    }

    /**
     * Required method for Iterator API.
     *
     * @throws UnsupportedOperationException
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove() not supported by FilteringIterator");
    }

    public void close() {
        CloserUtil.close(iterator);
    }

    /**
     * Gets the next record from the underlying iterator that passes the filter
     *
     * @return SAMRecord    the next filter-passing record
     */
    private SAMRecord getNextRecord() {
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (!filter.filterOut(record)) {
                return next;
            }
        }
        return null;
    }
}
