/*
 * Copyright (c) 2009 The Broad Institute
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

package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.sam.SAMReadValidator;
import org.broadinstitute.sting.utils.sam.SAMReadValidationException;
import org.broadinstitute.sting.utils.sam.SAMReadViolationHistogram;

import java.util.NoSuchElementException;

/**
 * A decorating iterator that examines the stream of reads, discarding those
 * that fail to meet a minimum standard for consumption by the GATK. 
 *
 * @author hanna
 * @version 0.1
 */

public class MalformedSAMFilteringIterator implements StingSAMIterator {
    /**
     * The wrapped iterator.  Get reads from here.
     */
    private StingSAMIterator wrapped = null;

    /**
     * Collector for SAM read violations.
     */
    private SAMReadViolationHistogram violations = null;

    /**
     * The next SAMRecord to return.;
     */
    private SAMRecord next = null;

    /**
     * Creates a new MalformedSAMFilteringIterator, and provides a collector for the count
     * @param wrapped The wrapped iterator to use as backing data.
     * @param violations A structure to hold a breakdown of validator violations.
     */
    public MalformedSAMFilteringIterator( StingSAMIterator wrapped, SAMReadViolationHistogram violations ) {
        this.wrapped = wrapped;
        this.violations = violations;
        seedNext();        
    }

    /**
     * Returns source information about the reads.
     * @return
     */
    public Reads getSourceInfo() {
        return wrapped.getSourceInfo();
    }

    /**
     * Gets an iterator, helpful for foreach loops.
     * @return An iterator sharing the same state variables as the current iterator.
     */
    public StingSAMIterator iterator() {
        return this;
    }

    /**
     * Checks to see whether there's a
     * @return True if a next is available, false otherwise.
     */
    public boolean hasNext() {
        return next != null;
    }

    /**
     * Gets the next valid record from the stream.
     * @return Next valid record.
     */
    public SAMRecord next() {
        SAMRecord current = next;
        if( current == null )
            throw new NoSuchElementException("MalformedSAMFilteringIterator: supply of reads is exhausted.");
        seedNext();
        return current;
    }

    /**
     * Closes the wrapped iterator.
     */
    public void close() {
        wrapped.close();
    }

    /**
     * Looks ahead for the next valid SAMRecord.
     */
    protected void seedNext() {
        next = null;        
        while( wrapped.hasNext() && next == null ) {
            SAMRecord toTest = wrapped.next();
            try {
                SAMReadValidator.validate(toTest);
                next = toTest;
            }
            catch ( SAMReadValidationException ex ) {
                violations.addViolation(ex);
            }
        }
    }

    /**
     * Throws an exception.  Remove is not supported.
     */
    public void remove() { throw new UnsupportedOperationException("Unable to remove from a StingSAMIterator"); }
}
