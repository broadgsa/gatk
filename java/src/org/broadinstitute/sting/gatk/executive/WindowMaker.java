package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.Reads;

import java.util.List;
import java.util.Queue;
import java.util.ArrayDeque;
import java.util.Iterator;

import net.sf.samtools.SAMRecord;
import net.sf.picard.util.PeekableIterator;

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
     * The data source for reads.  Will probably come directly from the BAM file.
     */
    private final StingSAMIterator sourceIterator;

    /**
     * Stores the sequence of intervals that the windowmaker should be tracking.
     */
    private final PeekableIterator<GenomeLoc> intervalIterator;

    /**
     * Which reads should be saved to go into the next interval?
     */
    private Queue<SAMRecord> overlappingReads = new ArrayDeque<SAMRecord>();

    /**
     * Create a new window maker with the given iterator as a data source, covering
     * the given inteervals.
     * @param iterator The data source for this window.
     * @param intervals The set of intervals over which to traverse.
     */
    public WindowMaker(StingSAMIterator iterator, List<GenomeLoc> intervals) {
        this.sourceIterator = iterator;
        this.intervalIterator = new PeekableIterator<GenomeLoc>(intervals.iterator());
    }

    public Iterator<WindowMakerIterator> iterator() {
        return this;
    }

    public boolean hasNext() {
        return intervalIterator.hasNext();
    }

    public WindowMakerIterator next() {
        return new WindowMakerIterator(intervalIterator.next());
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a window maker.");
    }

    public void close() {
        this.sourceIterator.close();
    }

    public class WindowMakerIterator implements StingSAMIterator {
        /**
         * The locus for which this iterator is currently returning reads.
         */
        private final GenomeLoc locus;

        /**
         * Which reads should be saved to go into the next interval?
         */
        private final Queue<SAMRecord> pendingOverlaps = new ArrayDeque<SAMRecord>();

        public WindowMakerIterator(GenomeLoc locus) {
            this.locus = locus;
        }

        public Reads getSourceInfo() {
            return sourceIterator.getSourceInfo();
        }

        public GenomeLoc getLocus() {
            return locus;
        }

        public WindowMakerIterator iterator() {
            return this;
        }

        public boolean hasNext() {
            return overlappingReads.size() > 0 || sourceIterator.hasNext();
        }

        public SAMRecord next() {
            SAMRecord nextRead = overlappingReads.size() > 0 ? overlappingReads.remove() : sourceIterator.next();
            if(intervalIterator.hasNext() && nextRead.getAlignmentEnd() >= intervalIterator.peek().getStart())
                pendingOverlaps.add(nextRead);
            return nextRead;
        }

        public void close() {
            overlappingReads = pendingOverlaps;    
        }

        public void remove() {
            throw new UnsupportedOperationException("Unable to remove from a window maker iterator.");
        }       
    }
}
