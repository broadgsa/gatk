/*
  * The Broad Institute
  * SOFTWARE COPYRIGHT NOTICE AGREEMENT
  * This software and its documentation are copyright Jan 22, 2009 by the
  * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
  *
  * This software is supplied without any warranty or guaranteed support whatsoever. Neither
  * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
  */
package edu.mit.broad.picard.sam;

import edu.mit.broad.picard.util.PeekableIterator;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;

import java.util.Comparator;

/**
 * Iterator for SAM records that implements comparable to enable sorting of iterators.
 * The comparison is performed by comparing the next record in the iterator to the next
 * record in another iterator and returning the ordering between those SAM records.
 */
class ComparableSamRecordIterator extends PeekableIterator<SAMRecord> implements Comparable<ComparableSamRecordIterator> {
    private Comparator<SAMRecord> comparator;
    private SAMFileReader reader;

    /**
     * Constructs an iterator for iteration over the supplied SAM file that will be
     * able to compare itself to other ComparableSAMRecordIterator instances using
     * the supplied comparator for ordering SAMRecords.
     *
     * @param sam the SAM file to read records from
     * @param comparator the Comparator to use to provide ordering fo SAMRecords
     */
    public ComparableSamRecordIterator(SAMFileReader sam, Comparator<SAMRecord> comparator) {
        super(sam.iterator());
        this.reader = sam;
        this.comparator = comparator;
    }

    /** Returns the reader from which this iterator was constructed. */
    public SAMFileReader getReader() {
        return reader;
    }

    /**
     * Compares this iterator to another comparable iterator based on the next record
     * available in each iterator.  If the two comparable iterators have different
     * comparator types internally an exception is thrown.
     *
     * @param that another iterator to compare to
     * @return a negative, 0 or positive number as described in the Comparator interface
     */
    public int compareTo(ComparableSamRecordIterator that) {
        if (this.comparator.getClass() != that.comparator.getClass()) {
            throw new IllegalStateException("Attempt to compare two ComparableSAMRecordIterators that " +
                    "have different orderings internally");
        }

        SAMRecord record = this.peek();
        SAMRecord record2 = that.peek();
        return comparator.compare(record, record2);
    }
}
