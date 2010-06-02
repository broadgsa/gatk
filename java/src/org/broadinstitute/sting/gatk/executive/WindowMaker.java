package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.filters.CountingFilteringIterator;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

import net.sf.samtools.SAMRecord;
import net.sf.picard.util.PeekableIterator;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;

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
    private final Reads sourceInfo;

    /**
     * Hold the read iterator so that it can be closed later.
     */
    private final StingSAMIterator readIterator;

    /**
     * The locus overflow tracker.
     */
    private final LocusOverflowTracker locusOverflowTracker;

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
     */
    public WindowMaker(StingSAMIterator iterator, List<GenomeLoc> intervals, List<SamRecordFilter> filters, List<LocusIteratorFilter> discards ) {
        this.sourceInfo = iterator.getSourceInfo();
        this.readIterator = iterator;

        LocusIterator locusIterator;
        Iterator<SAMRecord> wrappedIterator = TraversalEngine.addMandatoryFilteringIterators(iterator, filters);
        if(sourceInfo.getDownsamplingMethod() != null &&
          (sourceInfo.getDownsamplingMethod().type == DownsampleType.EXPERIMENTAL_BY_SAMPLE || sourceInfo.getDownsamplingMethod().type == DownsampleType.EXPERIMENTAL_NAIVE_DUPLICATE_ELIMINATOR)) {
            if ( discards.size() > 0 )
                throw new StingException("Experimental downsampling iterator doesn't support base discarding at this point; complain to Matt Hanna");
            locusIterator = new DownsamplingLocusIteratorByState(wrappedIterator,sourceInfo);
        } else
            locusIterator = new LocusIteratorByState(wrappedIterator,sourceInfo, discards);

        this.locusOverflowTracker = locusIterator.getLocusOverflowTracker();

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

        public Reads getSourceInfo() {
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
            // TODO: Move the monolithic sharding iterator so that we don't have to special case here.
            return sourceIterator.hasNext() && (locus == null || sourceIterator.peek().getLocation().overlapsP(locus));
        }

        public AlignmentContext next() {
            if(!hasNext()) throw new NoSuchElementException("WindowMakerIterator is out of elements for this interval.");
            return sourceIterator.next();
        }

        public LocusOverflowTracker getLocusOverflowTracker() {
            return locusOverflowTracker;
        }

        public void seedNextLocus() {
            // locus == null when doing monolithic sharding.
            // TODO: Move the monolithic sharding iterator so that we don't have to special case here.
            if(locus == null) return;

            while(sourceIterator.hasNext() && sourceIterator.peek().getLocation().isBefore(locus))
                sourceIterator.next();                
        }
    }
}
