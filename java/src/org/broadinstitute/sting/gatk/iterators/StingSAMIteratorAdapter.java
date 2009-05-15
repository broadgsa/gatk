package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.StingException;

/**
 *
 * User: aaron
 * Date: May 13, 2009
 * Time: 6:33:15 PM
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
 * @date May 13, 2009
 * <p/>
 * Class StingSAMIteratorAdapter
 * <p/>
 * This class adapts other SAMRecord iterators to the StingSAMIterator
 */
public class StingSAMIteratorAdapter {  

    public static StingSAMIterator adapt(Reads sourceInfo, Iterator<SAMRecord> iter) {
        return new PrivateStringSAMIterator(sourceInfo, iter);
    }

    public static StingSAMIterator adapt(Reads sourceInfo, CloseableIterator<SAMRecord> iter) {
        return new PrivateStringSAMCloseableIterator(sourceInfo, iter);
    }

}


/**
 * this class wraps iterators<SAMRecord> in a StingSAMIterator, which means just adding the
 * methods that implement the iterable<> interface and the close() method from CloseableIterator
 */
class PrivateStringSAMIterator implements StingSAMIterator {
    private Reads sourceInfo = null;
    private Iterator<SAMRecord> iter = null;

    PrivateStringSAMIterator(Reads sourceInfo, Iterator<SAMRecord> iter) {
        this.sourceInfo = sourceInfo;
        this.iter = iter;
    }

    public Reads getSourceInfo() {
        if( sourceInfo == null )
            throw new StingException("Unable to provide source info for the reads.  Please upgrade to the new data sharding framework.");
        return sourceInfo;
    }

    public void close() {
        // do nothing, we can't close the iterator anyway.
    }

    public boolean hasNext() {
        return iter.hasNext();
    }

    public SAMRecord next() {
        return iter.next();
    }

    public void remove() {
        throw new UnsupportedOperationException("StingSAMIterator's don't allow remove()ing");
    }

    public Iterator<SAMRecord> iterator() {
        return iter;
    }
}


/**
 * this class wraps closeable iterators<SAMRecord> in a StingSAMIterator, which means adding the
 * methods that implement the iterable<> interface.
 */
class PrivateStringSAMCloseableIterator implements StingSAMIterator {
    private Reads sourceInfo = null;
    private CloseableIterator<SAMRecord> iter = null;

    PrivateStringSAMCloseableIterator(Reads sourceInfo, CloseableIterator<SAMRecord> iter) {
        this.sourceInfo = sourceInfo;
        this.iter = iter;
    }

    public Reads getSourceInfo() {
        if( sourceInfo == null )
            throw new StingException("Unable to provide source info for the reads.  Please upgrade to the new data sharding framework.");
        return sourceInfo;
    }    

    public void close() {
        iter.close();
    }

    public boolean hasNext() {
        return iter.hasNext();
    }

    public SAMRecord next() {
        return iter.next();
    }

    public void remove() {
        throw new UnsupportedOperationException("StingSAMIterator's don't allow remove()ing");
    }

    public Iterator<SAMRecord> iterator() {
        return iter;
    }
}

