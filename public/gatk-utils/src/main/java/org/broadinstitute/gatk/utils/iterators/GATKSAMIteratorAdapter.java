/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.iterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

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
 * Class GATKSAMIteratorAdapter
 * <p/>
 * This class adapts other SAMRecord iterators to the GATKSAMIterator
 */
public class GATKSAMIteratorAdapter {

    public static GATKSAMIterator adapt(Iterator<SAMRecord> iter) {
        return new PrivateStringSAMIterator(iter);
    }

    public static GATKSAMIterator adapt(CloseableIterator<SAMRecord> iter) {
        return new PrivateStringSAMCloseableIterator(iter);
    }

}


/**
 * this class wraps iterators<SAMRecord> in a GATKSAMIterator, which means just adding the
 * methods that implement the iterable<> interface and the close() method from CloseableIterator
 */
class PrivateStringSAMIterator implements GATKSAMIterator {
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
        throw new UnsupportedOperationException("GATKSAMIterator's don't allow remove()ing");
    }

    public Iterator<SAMRecord> iterator() {
        return iter;
    }
}


/**
 * this class wraps closeable iterators<SAMRecord> in a GATKSAMIterator, which means adding the
 * methods that implement the iterable<> interface.
 */
class PrivateStringSAMCloseableIterator implements GATKSAMIterator {
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
        throw new UnsupportedOperationException("GATKSAMIterator's don't allow remove()ing");
    }

    public Iterator<SAMRecord> iterator() {
        return iter;
    }
}

