package org.broadinstitute.sting.gatk.iterators;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 16, 2009
 * Time: 6:08:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class SamQueryIterator implements Iterator<SAMRecord> {

    SAMFileReader reader = null;

    Iterator<GenomeLoc> locIter = null;
    CloseableIterator<SAMRecord> recordIter = null;

    public SamQueryIterator( SAMFileReader reader, ArrayList<GenomeLoc> locs ) {
        this.reader = reader;

        // Our internal contract for the class guarantees that locIter and recordIter are never null.
        // Initialize them and seed them with empty data as necessary.
        if(locs != null) {
            // The user requested a specific set of locations, set up the iterators accordly.
            locIter = locs.iterator();
            recordIter = new NullCloseableIterator<SAMRecord>();
        }
        else {
            // The user requested traversal of the entire SAM file.  Handle that here.
            // TODO: This would be better handled as a completely separate iterator.
            locIter = new ArrayList<GenomeLoc>().iterator();
            recordIter = reader.iterator();
        }

        bumpToNextSAMRecord();
    }

    public boolean hasNext() {
        bumpToNextSAMRecord();
        return recordIter.hasNext();
    }

    public SAMRecord next() {
        bumpToNextSAMRecord();
        return recordIter.next();
    }

    /**
     * Bump the loc iterator to the next spot with a read.
     *
     * For simplicity's sake, bumpToNextSAMRecord() expects locIter and recordIter to be non-null, and
     * guarantees that locIter and recordIter will be non-null after the bump.
     */
    private void bumpToNextSAMRecord() {
        // If there's a record still waiting in the current iterator, do nothing.
        if( recordIter.hasNext() ) {
            return;
        }

        // Otherwise, find the next record.
        recordIter.close();
        while( locIter.hasNext() ) {
            GenomeLoc currentLoc = locIter.next();
            recordIter = reader.queryOverlapping( currentLoc.getContig(),
                                                  (int)currentLoc.getStart(),
                                                  (int)currentLoc.getStop() );
            if( recordIter.hasNext() )
                break;
        }
    }
    
    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    private class NullCloseableIterator<T> implements CloseableIterator<T> {
        public boolean hasNext() { return false; }
        public T next() { throw new java.util.NoSuchElementException(); }
        public void close() {}
        public void remove() {}
    }
}
