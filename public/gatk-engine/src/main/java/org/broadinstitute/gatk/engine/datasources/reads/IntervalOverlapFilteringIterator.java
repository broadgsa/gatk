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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;

import java.util.List;
import java.util.NoSuchElementException;

/**
 * High efficiency filtering iterator designed to filter out reads only included
 * in the query results due to the granularity of the BAM index.
 *
 * Built into the BAM index is a notion of 16kbase granularity -- an index query for
 * two regions contained within a 16kbase chunk (say, chr1:5-10 and chr1:11-20) will
 * return exactly the same regions within the BAM file.  This iterator is optimized
 * to subtract out reads which do not at all overlap the interval list passed to the
 * constructor.
 *
 * Example:
 * interval list: chr20:6-10
 * Reads that would pass through the filter: chr20:6-10, chr20:1-15, chr20:1-7, chr20:8-15.
 * Reads that would be discarded by the filter: chr20:1-5, chr20:11-15.
 */
class IntervalOverlapFilteringIterator implements CloseableIterator<SAMRecord> {
    /**
     * The wrapped iterator.
     */
    private CloseableIterator<SAMRecord> iterator;

    /**
     * The next read, queued up and ready to go.
     */
    private SAMRecord nextRead;

    /**
     * Rather than using the straight genomic bounds, use filter out only mapped reads.
     */
    private boolean keepOnlyUnmappedReads;

    /**
     * Custom representation of interval bounds.
     * Makes it simpler to track current position.
     */
    private int[] intervalContigIndices;
    private int[] intervalStarts;
    private int[] intervalEnds;

    /**
     * Position within the interval list.
     */
    private int currentBound = 0;

    public IntervalOverlapFilteringIterator(CloseableIterator<SAMRecord> iterator, List<GenomeLoc> intervals) {
        this.iterator = iterator;

        // Look at the interval list to detect whether we should worry about unmapped reads.
        // If we find a mix of mapped/unmapped intervals, throw an exception.
        boolean foundMappedIntervals = false;
        for(GenomeLoc location: intervals) {
            if(! GenomeLoc.isUnmapped(location))
                foundMappedIntervals = true;
            keepOnlyUnmappedReads |= GenomeLoc.isUnmapped(location);
        }


        if(foundMappedIntervals) {
            if(keepOnlyUnmappedReads)
                throw new ReviewedGATKException("Tried to apply IntervalOverlapFilteringIterator to a mixed of mapped and unmapped intervals.  Please apply this filter to only mapped or only unmapped reads");
            this.intervalContigIndices = new int[intervals.size()];
            this.intervalStarts = new int[intervals.size()];
            this.intervalEnds = new int[intervals.size()];
            int i = 0;
            for(GenomeLoc interval: intervals) {
                intervalContigIndices[i] = interval.getContigIndex();
                intervalStarts[i] = interval.getStart();
                intervalEnds[i] = interval.getStop();
                i++;
            }
        }

        advance();
    }

    public boolean hasNext() {
        return nextRead != null;
    }

    public SAMRecord next() {
        if(nextRead == null)
            throw new NoSuchElementException("No more reads left in this iterator.");
        SAMRecord currentRead = nextRead;
        advance();
        return currentRead;
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from an IntervalOverlapFilteringIterator");
    }


    public void close() {
        iterator.close();
    }

    private void advance() {
        nextRead = null;

        if(!iterator.hasNext())
            return;

        SAMRecord candidateRead = iterator.next();
        while(nextRead == null && (keepOnlyUnmappedReads || currentBound < intervalStarts.length)) {
            if(!keepOnlyUnmappedReads) {
                // Mapped read filter; check against GenomeLoc-derived bounds.
                if(readEndsOnOrAfterStartingBound(candidateRead)) {
                    // This read ends after the current interval begins.
                    // Promising, but this read must be checked against the ending bound.
                    if(readStartsOnOrBeforeEndingBound(candidateRead)) {
                        // Yes, this read is within both bounds.  This must be our next read.
                        nextRead = candidateRead;
                        break;
                    }
                    else {
                        // Oops, we're past the end bound.  Increment the current bound and try again.
                        currentBound++;
                        continue;
                    }
                }
            }
            else {
                // Found a -L UNMAPPED read. NOTE: this is different than just being flagged as unmapped! We're done.
                if(AlignmentUtils.isReadGenomeLocUnmapped(candidateRead)) {
                    nextRead = candidateRead;
                    break;
                }
            }

            // No more reads available.  Stop the search.
            if(!iterator.hasNext())
                break;

            // No reasonable read found; advance the iterator.
            candidateRead = iterator.next();
        }
    }

    /**
     * Check whether the read lies after the start of the current bound.  If the read is unmapped but placed, its
     * end will be distorted, so rely only on the alignment start.
     * @param read The read to position-check.
     * @return True if the read starts after the current bounds.  False otherwise.
     */
    private boolean readEndsOnOrAfterStartingBound(final SAMRecord read) {
        return
                // Read ends on a later contig, or...
                read.getReferenceIndex() > intervalContigIndices[currentBound] ||
                        // Read ends of this contig...
                        (read.getReferenceIndex() == intervalContigIndices[currentBound] &&
                                // either after this location, or...
                                (read.getAlignmentEnd() >= intervalStarts[currentBound] ||
                                        // read is unmapped but positioned and alignment start is on or after this start point.
                                        (read.getReadUnmappedFlag() && read.getAlignmentStart() >= intervalStarts[currentBound])));
    }

    /**
     * Check whether the read lies before the end of the current bound.
     * @param read The read to position-check.
     * @return True if the read starts after the current bounds.  False otherwise.
     */
    private boolean readStartsOnOrBeforeEndingBound(final SAMRecord read) {
        return
                // Read starts on a prior contig, or...
                read.getReferenceIndex() < intervalContigIndices[currentBound] ||
                        // Read starts on this contig and the alignment start is registered before this end point.
                        (read.getReferenceIndex() == intervalContigIndices[currentBound] && read.getAlignmentStart() <= intervalEnds[currentBound]);
    }
}
