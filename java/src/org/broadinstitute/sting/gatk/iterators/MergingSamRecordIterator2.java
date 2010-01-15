/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.iterators;

import net.sf.picard.PicardException;
import net.sf.picard.sam.ReservedTagConstants;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

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
public class MergingSamRecordIterator2 implements CloseableIterator<SAMRecord>, Iterable<SAMRecord>, QueryIterator {
    protected PriorityQueue<ComparableSamRecordIterator> pq = null;
    protected final SamFileHeaderMerger samHeaderMerger;
    protected final SAMFileHeader.SortOrder sortOrder;
    protected static Logger logger = Logger.getLogger(MergingSamRecordIterator2.class);
    private SAMRecord mNextRecord;
    protected boolean initialized = false;
    protected final Reads reads;
    protected static boolean warnedUserAboutSortOrder = false; // so we only warn the user once

    /**
     * Constructs a new merging iterator with the same set of readers and sort order as
     * provided by the header merger parameter.
     *
     * @param headerMerger the header to merge
     * @param reads        the reads pile
     */
    public MergingSamRecordIterator2(final SamFileHeaderMerger headerMerger, Reads reads) {
        this.samHeaderMerger = headerMerger;
        this.reads = reads;
        this.sortOrder = headerMerger.getMergedHeader().getSortOrder();
        this.pq = new PriorityQueue<ComparableSamRecordIterator>(samHeaderMerger.getReaders().size());

    }

    /**
     * we hold off initializing because we don't know if they plan to query or just iterate.  This function
     * gets called if the hasNext() determines we haven't been query-ed yet, and we want to setup the iterator to start
     * at the beginning of the read pile
     */
    private void lazyInitialization() {
        if (initialized) {
            throw new UnsupportedOperationException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();
        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            checkSortOrder(reader);

            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(this.samHeaderMerger, reader, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    /**
     * verify the sort order
     *
     * @param reader the reader to check
     */
    private void checkSortOrder(SAMFileReader reader) {
        if (this.sortOrder != SAMFileHeader.SortOrder.unsorted && reader.getFileHeader().getSortOrder() != this.sortOrder) {
            String msg = String.format("The GATK requires your bam have %s sort order, but your BAM file header %s.  Continuing beyond this point is unsafe -- please update your BAM file to have a compatible sort order using samtools sort or Picard MergeBamFiles",
                    this.sortOrder, reader.getFileHeader().getAttribute("SO") == null ? "is missing the SO sort order flag" : "has an SO flag set to " + reader.getFileHeader().getAttribute("SO"));
            if (!reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.ALLOW_UNSET_BAM_SORT_ORDER)) {
                throw new PicardException(msg);
            } else if (!warnedUserAboutSortOrder) {
                warnedUserAboutSortOrder = true;
                Utils.warnUser(msg);
            }

        }
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
            checkSortOrder(reader);
            Iterator<SAMRecord> recordIter = reader.queryOverlapping(contig, start, stop);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(this.samHeaderMerger, reader, recordIter, comparator);
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
            checkSortOrder(reader);
            Iterator<SAMRecord> recordIter = reader.query(contig, start, stop, contained);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(this.samHeaderMerger, reader, recordIter, comparator);
            addIfNotEmpty(iterator);
        }
        setInitialized();

    }

    public void queryUnmappedReads() {
        if (initialized) {
            throw new IllegalStateException("You cannot double initialize a MergingSamRecordIterator2");
        }
        final SAMRecordComparator comparator = getComparator();
        for (final SAMFileReader reader : samHeaderMerger.getReaders()) {
            Iterator<SAMRecord> recordIter = null;
            if (reader.hasIndex()) {
                recordIter = reader.queryUnmapped();
            } else {
                // HACK: Supporting completely unmapped BAM files is easy.  Let's do a quick check to make sure
                //       these BAMs aren't partially indexed.
                if (reader.getFileHeader().getSequenceDictionary().size() > 0)
                    throw new StingException("Partially mapped BAM files without indices are not supported");
                recordIter = reader.iterator();
            }
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(this.samHeaderMerger, reader, recordIter, comparator);
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
            checkSortOrder(reader);
            Iterator<SAMRecord> recordIter = reader.queryContained(contig, start, stop);
            final ComparableSamRecordIterator iterator = new ComparableSamRecordIterator(this.samHeaderMerger, reader, recordIter, comparator);
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

    /** @return the next record, without moving the iterator */
    public SAMRecord peek() {
        return mNextRecord;
    }

    /**
     * Returns the next record from the top most iterator during merging.
     *
     * @return the next record, null if it's unavailable
     */
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

        // Fix the read group if needs be
        if (this.samHeaderMerger.hasReadGroupCollisions()) {
            final String oldGroupId = (String) record.getAttribute(ReservedTagConstants.READ_GROUP_ID);
            if (oldGroupId != null ) {
                final String newGroupId = this.samHeaderMerger.getReadGroupId(iterator.getReader(), oldGroupId);
                record.setAttribute(ReservedTagConstants.READ_GROUP_ID, newGroupId);
                }
        }

        // Fix the program group if needs be
        if (this.samHeaderMerger.hasProgramGroupCollisions()) {
            final String oldGroupId = (String) record.getAttribute(ReservedTagConstants.PROGRAM_GROUP_ID);
            if (oldGroupId != null ) {
                final String newGroupId = this.samHeaderMerger.getProgramGroupId(iterator.getReader(), oldGroupId);
                record.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID, newGroupId);
            }
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
     *
     * @param iterator the iterator to add
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
     *
     * @return the SAMRecordComparator appropriate for our sort order
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

    /**
     * Returns the merged header that the merging iterator is working from.
     *
     * @return our sam file header
     */
    public SAMFileHeader getHeader() {
        return this.samHeaderMerger.getMergedHeader();
    }


    /**
     * closes all open iterators....DO THIS or you will run out of handles
     * with sharding.
     */
    public void close() {
        for (ComparableSamRecordIterator iterator : pq)
            iterator.close();
    }


    /**
     * allows us to be used in the new style for loops
     *
     * @return this iterator, which can be used in for loop each iterations
     */
    public Iterator<SAMRecord> iterator() {
        return this;
    }

    /**
     * Gets source information for the reads.  Contains information about the original reads
     * files, plus information about downsampling, etc.
     *
     * @return the read object we were created for
     */
    public Reads getSourceInfo() {
        return this.reads;
    }
}

// Should replace picard class with the same name
class ComparableSamRecordIterator extends PeekableIterator<SAMRecord> implements Comparable<ComparableSamRecordIterator>, StingSAMIterator {
    private Reads sourceInfo;
    private final Comparator<SAMRecord> comparator;
    private final SAMFileReader reader;
    private final SamFileHeaderMerger mHeaderMerger;

    /**
     * Constructs an iterator for iteration over the supplied SAM file that will be
     * able to compare itself to other ComparableSAMRecordIterator instances using
     * the supplied comparator for ordering SAMRecords.
     *
     * @param sam        the SAM file to read records from
     * @param comparator the Comparator to use to provide ordering fo SAMRecords
     */
    public ComparableSamRecordIterator(SamFileHeaderMerger samHeaderMerger, final SAMFileReader sam, final Comparator<SAMRecord> comparator) {
        super(sam.iterator());
        this.reader = sam;
        this.comparator = comparator;
        mHeaderMerger = samHeaderMerger;
    }

    public ComparableSamRecordIterator(SamFileHeaderMerger samHeaderMerger, final SAMFileReader sam, Iterator<SAMRecord> iterator, final Comparator<SAMRecord> comparator) {
        super(iterator);    // use the provided iterator
        this.reader = sam;
        this.comparator = comparator;
        mHeaderMerger = samHeaderMerger;
    }

    public Reads getSourceInfo() {
        if (sourceInfo == null)
            throw new StingException("Unable to provide source info for the reads.  Please upgrade to the new data sharding framework.");
        return sourceInfo;
    }

    /**
     * Returns the reader from which this iterator was constructed.
     *
     * @return the SAMFileReader
     */
    public SAMFileReader getReader() {
        return reader;
    }

    /**
     * Compares this iterator to another comparable iterator based on the next record
     * available in each iterator.  If the two comparable iterators have different
     * comparator types internally an exception is thrown.
     *
     * @param that another iterator to compare to
     *
     * @return a negative, 0 or positive number as described in the Comparator interface
     */
    public int compareTo(final ComparableSamRecordIterator that) {
        if (this.comparator.getClass() != that.comparator.getClass()) {
            throw new IllegalStateException("Attempt to compare two ComparableSAMRecordIterators that " +
                    "have different orderings internally");
        }

        final SAMRecord record = this.peek();
        final SAMRecord record2 = that.peek();
        record.setHeader(mHeaderMerger.getMergedHeader());
        record2.setHeader(mHeaderMerger.getMergedHeader());
        int index, index2;
        try {
            index = mHeaderMerger.getMergedHeader().getSequenceIndex(record.getReferenceName());
            record.setReferenceIndex(index);

            index2 = mHeaderMerger.getMergedHeader().getSequenceIndex(record2.getReferenceName());
            record2.setReferenceIndex(index2);
        } catch (Exception e) {
            throw new StingException("MergingSamRecordIterator2: unable to correct the reference index for read " + record.getReadName() + " or record " + record2.getReadName(),e);
        }
        return comparator.compare(record, record2);
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
