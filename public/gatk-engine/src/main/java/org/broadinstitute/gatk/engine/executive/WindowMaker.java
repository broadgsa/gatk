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

package org.broadinstitute.gatk.engine.executive;

import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecordIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.locusiterator.LocusIterator;
import org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Transforms an iterator of reads which overlap the given interval list into an iterator of covered single-base loci
 * completely contained within the interval list.  To do this, it creates a LocusIteratorByState which will emit a single-bp
 * locus for every base covered by the read iterator, then uses the WindowMakerIterator.advance() to filter down that stream of
 * loci to only those covered by the given interval list.
 *
 * Example:
 * Incoming stream of reads: A:chr20:1-5, B:chr20:2-6, C:chr20:2-7, D:chr20:3-8, E:chr20:5-10
 * Incoming intervals: chr20:3-7
 *
 * Locus iterator by state will produce the following stream of data:
 *  chr1:1 {A}, chr1:2 {A,B,C}, chr1:3 {A,B,C,D}, chr1:4 {A,B,C,D}, chr1:5 {A,B,C,D,E},
 *  chr1:6 {B,C,D,E}, chr1:7 {C,D,E}, chr1:8 {D,E}, chr1:9 {E}, chr1:10 {E}
 *
 * WindowMakerIterator will then filter the incoming stream, emitting the following stream:
 *  chr1:3 {A,B,C,D}, chr1:4 {A,B,C,D}, chr1:5 {A,B,C,D,E}, chr1:6 {B,C,D,E}, chr1:7 {C,D,E}
 *
 * @author mhanna
 * @version 0.1
 */
public class WindowMaker implements Iterable<WindowMaker.WindowMakerIterator>, Iterator<WindowMaker.WindowMakerIterator> {
    /**
     * Source information for iteration.
     */
    private final ReadProperties sourceInfo;

    /**
     * Hold the read iterator so that it can be closed later.
     */
    private final GATKSAMRecordIterator readIterator;

    /**
     * The data source for reads.  Will probably come directly from the BAM file.
     */
    private final PeekableIterator<AlignmentContext> sourceIterator;

    /**
     * Stores the sequence of intervals that the windowmaker should be tracking.
     */
    private final PeekableIterator<GenomeLoc> intervalIterator;

    /**
     * In the case of monolithic sharding, this case returns whether the only shard has been generated.
     */
    private boolean shardGenerated = false;

    /**
     * The alignment context to return from this shard's iterator.  Lazy implementation: the iterator will not find the
     * currentAlignmentContext until absolutely required to do so.   If currentAlignmentContext is null and advance()
     * doesn't populate it, no more elements are available.  If currentAlignmentContext is non-null, currentAlignmentContext
     * should be returned by next().
     */
    private AlignmentContext currentAlignmentContext;

    /**
     * Create a new window maker with the given iterator as a data source, covering
     * the given intervals.
     * @param iterator The data source for this window.
     * @param intervals The set of intervals over which to traverse.
     * @param sampleNames The complete set of sample names in the reads in shard
     */

    private final LocusIteratorByState libs;

    public WindowMaker(Shard shard, GenomeLocParser genomeLocParser, GATKSAMIterator iterator, List<GenomeLoc> intervals, Collection<String> sampleNames) {
        this.sourceInfo = shard.getReadProperties();
        this.readIterator = new GATKSAMRecordIterator(iterator);

        this.libs = new LocusIteratorByState(readIterator,
                sourceInfo.getDownsamplingMethod(), sourceInfo.includeReadsWithDeletionAtLoci(),
                sourceInfo.keepUniqueReadListInLIBS(), genomeLocParser,sampleNames);
        this.sourceIterator = new PeekableIterator<AlignmentContext>(libs);

        this.intervalIterator = intervals.size()>0 ? new PeekableIterator<GenomeLoc>(intervals.iterator()) : null;
    }

    public WindowMaker(Shard shard, GenomeLocParser genomeLocParser, GATKSAMIterator iterator, List<GenomeLoc> intervals ) {
        this(shard, genomeLocParser, iterator, intervals, LocusIteratorByState.sampleListForSAMWithoutReadGroups());
    }

    public Iterator<WindowMakerIterator> iterator() {
        return this;
    }

    public boolean hasNext() {
        return (intervalIterator != null && intervalIterator.hasNext()) || !shardGenerated;
    }

    public WindowMakerIterator next() {
        shardGenerated = true;
        return new WindowMakerIterator(intervalIterator != null ? intervalIterator.next() : null);
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a window maker.");
    }

    public void close() {
        this.readIterator.close();
    }

    public class WindowMakerIterator extends LocusIterator {
        /**
         * The locus for which this iterator is currently returning reads.
         */
        private final GenomeLoc locus;

        public WindowMakerIterator(GenomeLoc locus) {
            this.locus = locus;
            advance();
        }

        public ReadProperties getSourceInfo() {
            return sourceInfo;
        }

        public GenomeLoc getLocus() {
            return locus;
        }

        public WindowMakerIterator iterator() {
            return this;
        }

        public boolean hasNext() {
            advance();
            return currentAlignmentContext != null;
        }

        public AlignmentContext next() {
            if(!hasNext()) throw new NoSuchElementException("WindowMakerIterator is out of elements for this interval.");

            // Consume this alignment context.
            AlignmentContext toReturn = currentAlignmentContext;
            currentAlignmentContext = null;

            // Return the current element.
            return toReturn;
        }

        private void advance() {
            // Need to find the next element that is not past shard boundaries.  If we travel past the edge of
            // shard boundaries, stop and let the next interval pick it up.
            while(currentAlignmentContext == null && sourceIterator.hasNext()) {
                // Advance the iterator and try again.
                AlignmentContext candidateAlignmentContext = sourceIterator.peek();

                if(locus == null) {
                    // No filter present.  Return everything that LocusIteratorByState provides us.
                    currentAlignmentContext = sourceIterator.next();
                }
                else if(locus.isPast(candidateAlignmentContext.getLocation()))
                    // Found a locus before the current window; claim this alignment context and throw it away.
                    sourceIterator.next();
                else if(locus.containsP(candidateAlignmentContext.getLocation())) {
                    // Found a locus within the current window; claim this alignment context and call it the next entry.
                    currentAlignmentContext = sourceIterator.next();
                }
                else if(locus.isBefore(candidateAlignmentContext.getLocation())) {
                    // Whoops.  Skipped passed the end of the region.  Iteration for this window is complete.  Do
                    // not claim this alignment context in case it is part of the next shard.
                    break;
                }
                else
                    throw new ReviewedGATKException("BUG: filtering locus does not contain, is not before, and is not past the given alignment context");
            }
        }

        @Override
        public LocusIteratorByState getLIBS() {
            return libs;
        }
    }
}
