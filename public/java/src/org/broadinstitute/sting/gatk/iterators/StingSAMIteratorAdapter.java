package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

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

    public static StingSAMIterator adapt(Iterator<SAMRecord> iter) {
        return new PrivateStringSAMIterator(iter);
    }

    public static StingSAMIterator adapt(CloseableIterator<SAMRecord> iter) {
        return new PrivateStringSAMCloseableIterator(iter);
    }

}


/**
 * this class wraps iterators<SAMRecord> in a StingSAMIterator, which means just adding the
 * methods that implement the iterable<> interface and the close() method from CloseableIterator
 */
class PrivateStringSAMIterator implements StingSAMIterator {
    private Iterator<SAMRecord> iter = null;

    PrivateStringSAMIterator(Iterator<SAMRecord> iter) {
        this.iter = iter;
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
    private CloseableIterator<SAMRecord> iter = null;

    PrivateStringSAMCloseableIterator(CloseableIterator<SAMRecord> iter) {
        this.iter = iter;
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

