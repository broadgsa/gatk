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
