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
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.util.Iterator;


public class LegacyDownsampleIterator implements StingSAMIterator {

    StingSAMIterator it;
    int cutoff;
    SAMRecord next;

    public LegacyDownsampleIterator(StingSAMIterator it, double fraction) {
        this.it = it;
        cutoff = (int)(fraction * 10000);
        next = getNextRecord();
    }

    public boolean hasNext() {
        return next != null;
    }

    public SAMRecord next()  {
        SAMRecord result = next;
        next = getNextRecord();
        return result;
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    private SAMRecord getNextRecord() {
        while ( true ) {
            if ( !it.hasNext() )
                return null;
            SAMRecord rec = it.next();
            if ( GenomeAnalysisEngine.getRandomGenerator().nextInt(10000) < cutoff )
                return rec;
        }
    }

    public void close() {
        it.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}