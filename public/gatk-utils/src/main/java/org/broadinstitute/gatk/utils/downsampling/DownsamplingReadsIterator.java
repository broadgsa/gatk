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

package org.broadinstitute.gatk.utils.downsampling;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;


/**
 * GATKSAMIterator wrapper around our generic reads downsampler interface. Converts the push-style
 * downsampler interface to a pull model.
 *
 * @author David Roazen
 */
public class DownsamplingReadsIterator implements GATKSAMIterator {

    private GATKSAMIterator nestedSAMIterator;
    private ReadsDownsampler<SAMRecord> downsampler;
    private Collection<SAMRecord> downsampledReadsCache;
    private SAMRecord nextRead = null;
    private Iterator<SAMRecord> downsampledReadsCacheIterator = null;

    /**
     * @param iter wrapped iterator from which this iterator will pull reads
     * @param downsampler downsampler through which the reads will be fed
     */
    public DownsamplingReadsIterator( GATKSAMIterator iter, ReadsDownsampler<SAMRecord> downsampler ) {
        nestedSAMIterator = iter;
        this.downsampler = downsampler;

        advanceToNextRead();
    }

    public boolean hasNext() {
        return nextRead != null;
    }

    public SAMRecord next() {
        if ( nextRead == null ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        SAMRecord toReturn = nextRead;
        advanceToNextRead();

        return toReturn;
    }

    private void advanceToNextRead() {
        if ( ! readyToReleaseReads() && ! fillDownsampledReadsCache() ) {
            nextRead = null;
        }
        else {
            nextRead = downsampledReadsCacheIterator.next();
        }
    }

    private boolean readyToReleaseReads() {
        return downsampledReadsCacheIterator != null && downsampledReadsCacheIterator.hasNext();
    }

    private boolean fillDownsampledReadsCache() {
        while ( nestedSAMIterator.hasNext() && ! downsampler.hasFinalizedItems() ) {
            downsampler.submit(nestedSAMIterator.next());
        }

        if ( ! nestedSAMIterator.hasNext() ) {
            downsampler.signalEndOfInput();
        }

        // use returned collection directly rather than make a copy, for speed
        downsampledReadsCache = downsampler.consumeFinalizedItems();
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