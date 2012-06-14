/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;


/**
 * StingSAMIterator wrapper around our generic reads downsampler interface
 *
 * @author David Roazen
 */
public class DownsamplingReadsIterator implements StingSAMIterator {

    private StingSAMIterator nestedSAMIterator;
    private ReadsDownsampler<SAMRecord> downsampler;
    private Collection<SAMRecord> downsampledReadsCache;
    private Iterator<SAMRecord> downsampledReadsCacheIterator;

    public DownsamplingReadsIterator( StingSAMIterator iter, ReadsDownsampler<SAMRecord> downsampler ) {
        nestedSAMIterator = iter;
        this.downsampler = downsampler;
        fillDownsampledReadsCache();
    }

    public boolean hasNext() {
        if ( downsampledReadsCacheIterator.hasNext() ) {
            return true;
        }
        else if ( ! nestedSAMIterator.hasNext() || ! fillDownsampledReadsCache() ) {
            return false;
        }

        return true;
    }

    public SAMRecord next() {
        if ( ! downsampledReadsCacheIterator.hasNext() && ! fillDownsampledReadsCache() ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        return downsampledReadsCacheIterator.next();
    }

    private boolean fillDownsampledReadsCache() {
        while ( nestedSAMIterator.hasNext() && ! downsampler.hasDownsampledItems() ) {
            downsampler.submit(nestedSAMIterator.next());
        }

        if ( ! nestedSAMIterator.hasNext() ) {
            downsampler.signalEndOfInput();
        }

        downsampledReadsCache = downsampler.consumeDownsampledItems();
        downsampledReadsCacheIterator = downsampledReadsCache.iterator();

        return downsampledReadsCacheIterator.hasNext();
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() {
        nestedSAMIterator.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}