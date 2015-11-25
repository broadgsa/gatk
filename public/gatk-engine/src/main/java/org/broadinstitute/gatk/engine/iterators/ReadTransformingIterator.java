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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Iterator;

/**
 * Iterator that applies a ReadTransformer to a stream of reads
 */
public class ReadTransformingIterator implements GATKSAMIterator {
    private final GATKSAMIterator it;
    private final ReadTransformer transformer;

    /**
     * Creates a new ReadTransforming iterator
     */
    @Requires({"it != null", "transformer != null", "transformer.isInitialized()"})
    public ReadTransformingIterator(final GATKSAMIterator it, final ReadTransformer transformer) {
        if ( ! transformer.isInitialized() )
            throw new IllegalStateException("Creating a read transformer stream for an uninitialized read transformer: " + transformer);
        if ( transformer.getApplicationTime() == ReadTransformer.ApplicationTime.FORBIDDEN )
            throw new IllegalStateException("Creating a read transformer stream for a forbidden transformer " + transformer);

        this.it = it;
        this.transformer = transformer;
    }

    @Requires("hasNext()")
    @Ensures("result != null")
    public SAMRecord next()     {
        final GATKSAMRecord read = (GATKSAMRecord)it.next();
        return transformer.apply(read);
    }

    public boolean hasNext()    { return this.it.hasNext(); }
    public void remove()        { throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!"); }
    public void close()         { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
