/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Iterator;

/**
 * Temporarily hack to convert SAMRecords to GATKSAMRecords
 *
 * User: depristo
 * Date: 1/11/13
 * Time: 1:19 PM
 */
public class GATKSAMIterator implements CloseableIterator<GATKSAMRecord>, Iterable<GATKSAMRecord> {
    final CloseableIterator<SAMRecord> it;

    public GATKSAMIterator(final CloseableIterator<SAMRecord> it) {
        this.it = it;
    }

    public GATKSAMIterator(final StingSAMIterator it) {
        this.it = it;
    }

    @Override public boolean hasNext() { return it.hasNext(); }
    @Override public GATKSAMRecord next() { return (GATKSAMRecord)it.next(); }
    @Override public void remove() { it.remove(); }
    @Override public void close() { it.close(); }
    @Override public Iterator<GATKSAMRecord> iterator() { return this; }
}
