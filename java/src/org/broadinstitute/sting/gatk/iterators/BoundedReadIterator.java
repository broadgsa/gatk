package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import java.util.Iterator;

/**
 *
 * User: aaron
 * Date: Apr 14, 2009
 * Time: 5:27:31 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BoundedReadIterator
 * <p/>
 * This class implements a read iterator that is bounded by the number of reads
 * it will produce over the iteration.
 */
public class BoundedReadIterator implements StingSAMIterator {

    // the genome loc we're bounding
    final private long readCount;
    private long currentCount = 0;

    // the iterator we want to decorate
    private final StingSAMIterator iterator;

    // our unmapped read flag
    private boolean doNotUseThatUnmappedReadPile = false;

    // are we open
    private boolean isOpen = false;

    // the next read we've buffered
    private SAMRecord record = null;

    /**
     * constructor
     * @param iter
     * @param readCount
     */
    public BoundedReadIterator(StingSAMIterator iter, long readCount) {
        if (iter != null) {
            isOpen = true;

        }
        this.iterator = iter;
        this.readCount = readCount;
    }

    public void useUnmappedReads(boolean useThem) {
        this.doNotUseThatUnmappedReadPile = useThem;
    }


    public SAMFileHeader getHeader() {
        // todo: this is bad, we need an iterface out there for samrecords that supports getting the header,
        // regardless of the merging
        if (iterator instanceof MergingSamRecordIterator2)
            return ((MergingSamRecordIterator2)iterator).getHeader();
        else
            return null;
    }

    /**
     * Do we have a next? If the iterator has a read and we're not over the read
     * count, then yes
     * @return
     */
    public boolean hasNext() {
        if (isOpen && iterator.hasNext() && currentCount < readCount) {
            record = iterator.next();
            ++currentCount;
            if (record.getAlignmentStart() == 0 && doNotUseThatUnmappedReadPile) {
                return false;
            }
            return true;
        } else {
            if (isOpen) {
                close();
            }
            return false;
        }
    }

    /**
     * get the next SAMRecord
     * @return SAMRecord representing the next read
     */
    public SAMRecord next() {
        return record;
    }

    /**
     * this is unsupported on SAMRecord iterators
     */
    public void remove() {
        throw new UnsupportedOperationException("You cannot use an iterator to remove a SAMRecord");
    }

    /**
     * close the iterator
     */
    public void close() {
        isOpen = false;
        iterator.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
