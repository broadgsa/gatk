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
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.util.Iterator;

/**
 * Traps BAM formatting errors in underlying iterator and rethrows meaningful GATK UserExceptions
 */
public class MalformedBAMErrorReformatingIterator implements CloseableIterator<SAMRecord> {
    File source;
    CloseableIterator<SAMRecord> it;

    public MalformedBAMErrorReformatingIterator(final File source, final CloseableIterator<SAMRecord> it) {
        this.it = it;
        this.source = source;
    }

    public boolean hasNext() {
        try {
            return this.it.hasNext();
        } catch ( RuntimeException e ) { // we need to catch RuntimeExceptions here because the Picard code is throwing them (among SAMFormatExceptions) sometimes
            throw new UserException.MalformedBAM(source, e.getMessage());
        }
    }

    public SAMRecord next() {
        try {
            return it.next();
        } catch ( RuntimeException e ) { // we need to catch RuntimeExceptions here because the Picard code is throwing them (among SAMFormatExceptions) sometimes
            throw new UserException.MalformedBAM(source, e.getMessage());
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
