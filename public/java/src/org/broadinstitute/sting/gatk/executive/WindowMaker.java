package org.broadinstitute.sting.gatk.executive;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorByState;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Buffer shards of data which may or may not contain multiple loci into
 * iterators of all data which cover an interval.  Its existence is an homage
 * to Mark's stillborn WindowMaker, RIP 2009.
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
    private final StingSAMIterator readIterator;

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

    public WindowMaker(Shard shard, GenomeLocParser genomeLocParser, StingSAMIterator iterator, List<GenomeLoc> intervals, Collection<String> sampleNames) {
        this.sourceInfo = shard.getReadProperties();
        this.readIterator = iterator;
        this.sourceIterator = new PeekableIterator<AlignmentContext>(new LocusIteratorByState(iterator,sourceInfo,genomeLocParser, sampleNames));
        this.intervalIterator = intervals.size()>0 ? new PeekableIterator<GenomeLoc>(intervals.iterator()) : null;
    }

    public WindowMaker(Shard shard, GenomeLocParser genomeLocParser, StingSAMIterator iterator, List<GenomeLoc> intervals ) {
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
                    throw new ReviewedStingException("BUG: filtering locus does not contain, is not before, and is not past the given alignment context");
            }
        }
    }
}
