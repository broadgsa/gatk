/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever.
* Neither the Broad Institute nor MIT can be responsible for its use, misuse, or
* functionality.
*/
package org.broadinstitute.sting.gatk.iterators;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.sam.ReservedTagConstants;
import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import edu.mit.broad.picard.util.PeekableIterator;
import net.sf.samtools.*;
import org.apache.log4j.Logger;

import java.lang.reflect.Constructor;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Provides an iterator interface for merging multiple underlying iterators into a single
 * iterable stream. The underlying iterators/files must all have the same sort order unless
 * the requested output format is unsorted, in which case any combination is valid.
 */
public class MergingSamRecordIterator2 implements StingSAMIterator {
    protected PriorityQueue<ComparableSamRecordIterator> pq = null;
    protected final SamFileHeaderMerger samHeaderMerger;
    protected final SAMFileHeader.SortOrder sortOrder;
    protected static Logger logger = Logger.getLogger(MergingSamRecordIterator2.class);
    private SAMRecord mNextRecord;
    protected boolean initialized = false;

    /**
     * Constructs a new merging iterator with the same set of readers and sort order as
     * provided by the header merger parameter.
     */
    public MergingSamRecordIterator2(final SamFileHeaderMerger headerMerger) {
        this.samHeaderMerger = headerMerger;
        this.sortOrder = headerMerger.getMergedHeader().getSortOrder();
        this.pq = new PriorityQueue<ComparableSamRecordIterator>(samHeaderMerger.getReaders().size());

    }

    /** this class MUST only be initialized once, since the creation of the */
    private void lazyInitialization() {
        if (initialized) {
            throw new UnsupportedOperationException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();
        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            if (this.sortOrder != SAMFileHeader.SortOrder.unsorted && reader.getFileHeader().getSortOrder() != this.sortOrder) {
                throw new PicardException("Files are not compatible with sort order: " + reader.getFileHeader().getSortOrder() + " vrs " + this.sortOrder);
            }

            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(reader, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    public boolean supportsSeeking() {
        return true;
    }

    public void queryOverlapping(final String contig, final int start, final int stop) {
        if (initialized) {
            throw new IllegalStateException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();

        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            Iterator<SAMRecord> recordIter = reader.queryOverlapping(contig, start, stop);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(reader, recordIter, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    public void query(final String contig, final int start, final int stop, final boolean contained) {
        if (initialized) {
            throw new IllegalStateException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();
        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            //reader.close();
            Iterator<SAMRecord> recordIter = reader.query(contig, start, stop, contained);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(reader, recordIter, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    public void queryContained(final String contig, final int start, final int stop) {
        if (initialized) {
            throw new IllegalStateException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();
        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            Iterator<SAMRecord> recordIter = reader.queryContained(contig, start, stop);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(reader, recordIter, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    private void setInitialized() {
        initialized = true;
        mNextRecord = nextRecord();
    }


    /** Returns true if any of the underlying iterators has more records, otherwise false. */
    public boolean hasNext() {
        if (!initialized) {
            lazyInitialization();
        }
        if (this.pq.isEmpty() && mNextRecord == null) {
            return false;
        }
        return true;
    }

    public SAMRecord next() {
        SAMRecord r = mNextRecord;
        if (!this.pq.isEmpty()) {
            mNextRecord = nextRecord();
        } else {
            mNextRecord = null;
        }
        return r;
    }

    public SAMRecord peek() {
        return mNextRecord;
    }

    /** Returns the next record from the top most iterator during merging. */
    public SAMRecord nextRecord() {
        if (!initialized) {
            lazyInitialization();
        }

        final ComparableSamRecordIterator iterator = this.pq.poll();

        if (iterator == null) {
            return null;
        }
        final SAMRecord record = iterator.next();
        addIfNotEmpty(iterator);

        if (this.samHeaderMerger.hasGroupIdDuplicates()) {
            final String id = (String) record.getAttribute(ReservedTagConstants.READ_GROUP_ID);
            final String newId = this.samHeaderMerger.getReadGroupId(iterator.getReader(), id);
            record.setAttribute(ReservedTagConstants.READ_GROUP_ID, newId);
        }
        final String oldProgramGroupId = (String) record.getAttribute(SAMTag.PG.toString());
        if (oldProgramGroupId != null) {
            final String newProgramGroupId = this.samHeaderMerger.getProgramGroupId(iterator.getReader(), oldProgramGroupId);
            record.setAttribute(SAMTag.PG.toString(), newProgramGroupId);
        }

        // if we don't have a read group, set one correctly
        if (record.getAttribute(SAMTag.RG.toString()) == null) {
            List<SAMReadGroupRecord> readGroups = iterator.getReader().getFileHeader().getReadGroups();
            if (readGroups.size() == 1) {
                record.setAttribute(SAMTag.RG.toString(), readGroups.get(0).getReadGroupId());
                record.setAttribute(SAMTag.SM.toString(), readGroups.get(0).getReadGroupId());
            } else {
                logger.warn("Unable to set read group of ungrouped read: unable to pick default group, there are " + readGroups.size() + " possible.");
            }
        }

        record.setHeader(samHeaderMerger.getMergedHeader());
        //System.out.printf("NEXT = %s %s %d%n", record.getReadName(), record.getReferenceName(), record.getAlignmentStart());
        //System.out.printf("PEEK = %s %s %d%n", this.pq.peek().peek().getReadName(), this.pq.peek().peek().getReferenceName(), this.pq.peek().peek().getAlignmentStart());

        return record;
    }

    /**
     * Adds iterator to priority queue. If the iterator has more records it is added
     * otherwise it is closed and not added.
     */
    protected void addIfNotEmpty(final ComparableSamRecordIterator iterator) {
        //System.out.printf("Adding %s %s %d%n", iterator.peek().getReadName(), iterator.peek().getReferenceName(), iterator.peek().getAlignmentStart());
        if (iterator.hasNext()) {
            pq.offer(iterator);
        } else {
            iterator.close();
        }
    }

    /** Unsupported operation. */
    public void remove() {
        throw new UnsupportedOperationException("MergingSAMRecorderIterator.remove()");
    }

    /**
     * Get the right comparator for a given sort order (coordinate, alphabetic). In the
     * case of "unsorted" it will return a comparator that gives an arbitrary but reflexive
     * ordering.
     */
    protected SAMRecordComparator getComparator() {
        // For unsorted build a fake comparator that compares based on object ID
        if (this.sortOrder == SAMFileHeader.SortOrder.unsorted) {
            return new SAMRecordComparator() {
                public int fileOrderCompare(final SAMRecord lhs, final SAMRecord rhs) {
                    return System.identityHashCode(lhs) - System.identityHashCode(rhs);
                }

                public int compare(final SAMRecord lhs, final SAMRecord rhs) {
                    return fileOrderCompare(lhs, rhs);
                }
            };
        }

        // Otherwise try and figure out what kind of comparator to return and build it
        final Class<? extends SAMRecordComparator> type = this.sortOrder.getComparator();

        try {
            final Constructor<? extends SAMRecordComparator> ctor = type.getConstructor(SAMFileHeader.class);
            //System.out.printf("Getting comparator %s%n", ctor.toGenericString());
            return ctor.newInstance(this.samHeaderMerger.getMergedHeader());
        }
        catch (Exception e) {
            try {
                final Constructor<? extends SAMRecordComparator> ctor = type.getConstructor();
                return ctor.newInstance();
            }
            catch (Exception e2) {
                throw new PicardException("Could not instantiate a comparator for sort order: " + this.sortOrder, e2);
            }
        }
    }

    /** Returns the merged header that the merging iterator is working from. */
    public SAMFileHeader getHeader() {
        return this.samHeaderMerger.getMergedHeader();
    }


    /**
     * closes all the file handles for the readers....DO THIS or you will run out of handles
     * with sharding.
     */
    public void close() {
        for (SAMFileReader reader : samHeaderMerger.getReaders()) {
            reader.close();
        }
    }


    /**
     * allows us to be used in the new style for loops
     *
     * @return
     */
    public Iterator<SAMRecord> iterator() {
        return this;
    }
}

// Should replace picard class with the same name
class ComparableSamRecordIterator extends PeekableIterator<SAMRecord> implements Comparable<ComparableSamRecordIterator>, StingSAMIterator {
    private final Comparator<SAMRecord> comparator;
    private final SAMFileReader reader;

    /**
     * Constructs an iterator for iteration over the supplied SAM file that will be
     * able to compare itself to other ComparableSAMRecordIterator instances using
     * the supplied comparator for ordering SAMRecords.
     *
     * @param sam        the SAM file to read records from
     * @param comparator the Comparator to use to provide ordering fo SAMRecords
     */
    public ComparableSamRecordIterator(final SAMFileReader sam, final Comparator<SAMRecord> comparator) {
        super(sam.iterator());
        this.reader = sam;
        this.comparator = comparator;
    }

    public ComparableSamRecordIterator(final SAMFileReader sam, Iterator<SAMRecord> iterator, final Comparator<SAMRecord> comparator) {
        super(iterator);    // use the provided iterator
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
    public int compareTo(final ComparableSamRecordIterator that) {
        if (this.comparator.getClass() != that.comparator.getClass()) {
            throw new IllegalStateException("Attempt to compare two ComparableSAMRecordIterators that " +
                    "have different orderings internally");
        }

        final SAMRecord record = this.peek();
        final SAMRecord record2 = that.peek();
        //System.out.printf("Comparing %s vs. %s => %d%n", record.getReadName(), record2.getReadName(), comparator.compare(record, record2));
        return comparator.compare(record, record2);
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
