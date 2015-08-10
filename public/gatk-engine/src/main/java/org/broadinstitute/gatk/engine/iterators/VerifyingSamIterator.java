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

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

import java.util.Iterator;

/**
 * Verifies that the incoming stream of reads is correctly sorted
 */
public class VerifyingSamIterator implements GATKSAMIterator {
    GATKSAMIterator it;
    SAMRecord last = null;
    boolean checkOrderP = true;

    public VerifyingSamIterator(GATKSAMIterator it) {
        this.it = it;
    }

    public boolean hasNext() { return this.it.hasNext(); }
    public SAMRecord next() {

        SAMRecord cur = it.next();
        if ( last != null )
            verifyRecord(last, cur);
        if ( ! cur.getReadUnmappedFlag() )
            last = cur;
        return cur;
    }

    private void verifyRecord( final SAMRecord last, final SAMRecord cur ) {
        if ( checkOrderP && isOutOfOrder(last, cur) ) {
            this.last = null;
            throw new UserException.MissortedBAM(String.format("reads are out of order:%nlast:%n%s%ncurrent:%n%s%n", last.format(), cur.format()) );
        }
    }

    private boolean isOutOfOrder( final SAMRecord last, final SAMRecord cur ) {
        if ( last == null || cur.getReadUnmappedFlag() )
            return false;
        else {
            if(last.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX || last.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START)
                throw new UserException.MalformedBAM(last,String.format("read %s has inconsistent mapping information.",last.format()));
            if(cur.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX || cur.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START)
                throw new UserException.MalformedBAM(last,String.format("read %s has inconsistent mapping information.",cur.format()));

            return (last.getReferenceIndex() > cur.getReferenceIndex()) ||
                    (last.getReferenceIndex().equals(cur.getReferenceIndex()) &&
                            last.getAlignmentStart() > cur.getAlignmentStart());
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() {
        it.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
