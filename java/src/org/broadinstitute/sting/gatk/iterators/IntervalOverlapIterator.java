package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Iterator;


/**
 * 
 * @author aaron 
 * 
 * Class DuplicateDetectorIterator
 *
 * remove reads that overlap the passed in interval, yea!
 *
 */
public class IntervalOverlapIterator implements StingSAMIterator {

    // our wrapped iterator
    private final StingSAMIterator mIter;
    private final boolean throwException;

    // storage for the next record
    private SAMRecord mNextRecord = null;

    // the genomic location which we filter on
    private final GenomeLoc mLoc;

    /**
     * Create a DuplicateDetectorIterator from another sam iterator
     * @param iter something that implements StingSAMIterator
     * @param blowUpOnDup if we find a dup, do we throw an exception (blow up) or do we drop it
     */
    public IntervalOverlapIterator(StingSAMIterator iter, GenomeLoc filterLocation, boolean blowUpOnDup) {
        this.mIter = iter;
        this.throwException = blowUpOnDup;
        this.mLoc = filterLocation;
        if (iter.hasNext()) {
            next();
        }
    }

    /**
     * Gets source information for the reads.  Contains information about the original reads
     * files, plus information about downsampling, etc.
     *
     * @return
     */
    public Reads getSourceInfo() {
        return mIter.getSourceInfo();
    }

    /**
     * close this iterator
     */
    public void close() {
        if (mIter != null) mIter.close();
    }

    /**
     * do we have a next?
     * @return true if yes, false if not
     */
    public boolean hasNext() {
        return (mNextRecord != null);
    }

    /**
     * get the next record
     * @return a SAMRecord
     */
    public SAMRecord next() {
        SAMRecord ret = mNextRecord;
        while (mIter.hasNext()) {
            mNextRecord = mIter.next();
            if (!isOverlaping(mNextRecord)) return ret;
        }
        mNextRecord = null;
        return ret;
    }

    /**
     * not supported
     */
    public void remove() {
        throw new UnsupportedOperationException("You can't call remove, so, like, I guess please don't");
    }

    /**
     * create an iterator out of the this type
     * @return this!
     */
    public Iterator<SAMRecord> iterator() {
        return this;
    }

    /**
     * determine if a read overlaps the specified interval that was passed in
     * @param rec the read
     * @return true if it overlaps, false otherwise
     */
    private boolean isOverlaping(SAMRecord rec) {
        return mLoc.overlapsP(GenomeLocParser.createGenomeLoc(rec));
    }

}
