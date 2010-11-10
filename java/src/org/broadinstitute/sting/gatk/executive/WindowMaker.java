package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.utils.GenomeLocParser;

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
     * Create a new window maker with the given iterator as a data source, covering
     * the given intervals.
     * @param iterator The data source for this window.
     * @param intervals The set of intervals over which to traverse.
     * @param discards a filter at that indicates read position relative to some locus?
     */
    public WindowMaker(Shard shard, GenomeLocParser genomeLocParser, StingSAMIterator iterator, List<GenomeLoc> intervals, List<LocusIteratorFilter> discards ) {
        this.sourceInfo = shard.getReadProperties();
        this.readIterator = iterator;

        LocusIterator locusIterator = new LocusIteratorByState(iterator,sourceInfo,genomeLocParser,discards);

        this.sourceIterator = new PeekableIterator<AlignmentContext>(locusIterator);
        this.intervalIterator = intervals.size()>0 ? new PeekableIterator<GenomeLoc>(intervals.iterator()) : null;
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
            seedNextLocus();
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
            // locus == null when doing monolithic sharding.
            return sourceIterator.hasNext() && sourceIterator.peek().getLocation().overlapsP(locus);
        }

        public AlignmentContext next() {
            if(!hasNext()) throw new NoSuchElementException("WindowMakerIterator is out of elements for this interval.");
            return sourceIterator.next();
        }

        public void seedNextLocus() {
            // locus == null when doing monolithic sharding.
            while(sourceIterator.hasNext() && sourceIterator.peek().getLocation().isBefore(locus))
                sourceIterator.next();                
        }
    }
}
