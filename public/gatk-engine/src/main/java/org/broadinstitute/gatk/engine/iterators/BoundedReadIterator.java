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

package org.broadinstitute.gatk.engine.iterators;

import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

import java.util.Iterator;

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
public class BoundedReadIterator implements GATKSAMIterator {

    // the genome loc we're bounding
    final private long readCount;
    private long currentCount = 0;

    // the iterator we want to decorate
    private final GATKSAMIterator iterator;

    // our unmapped read flag
    private boolean doNotUseThatUnmappedReadPile = false;

    /**
     * The next read that we've buffered.  Null indicates that there's
     * nothing in the buffer (not that there isn't a next read).
     */
    private SAMRecord record = null;

    /**
     * constructor
     * @param iter
     * @param readCount
     */
    public BoundedReadIterator(GATKSAMIterator iter, long readCount) {
        this.iterator = iter;
        this.readCount = readCount;
    }

    public void useUnmappedReads(boolean useThem) {
        this.doNotUseThatUnmappedReadPile = useThem;
    }

    public SAMFileHeader getHeader() {
        // todo: this is bad, we need an iterface out there for samrecords that supports getting the header,
        // regardless of the merging
        if (iterator instanceof MergingSamRecordIterator)
            return ((MergingSamRecordIterator)iterator).getMergedHeader();
        else
            return null;
    }

    /**
     * Do we have a next? If the iterator has a read and we're not over the read
     * count, then yes
     * @return
     */
    public boolean hasNext() {
        if( record != null )
            return true;

        if (iterator.hasNext() && currentCount < readCount) {
            record = iterator.next();
            ++currentCount;
            if (record.getAlignmentStart() == 0 && doNotUseThatUnmappedReadPile) {
                return false;
            }
            return true;
        } else {
            return false;
        }
    }

    /**
     * get the next SAMRecord
     * @return SAMRecord representing the next read
     */
    public SAMRecord next() {
        SAMRecord cached = record;
        record = null;
        return cached;
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
        iterator.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
