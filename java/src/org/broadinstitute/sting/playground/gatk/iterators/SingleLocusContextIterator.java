package org.broadinstitute.sting.playground.gatk.iterators;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.iterators.LocusContextIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Predicate;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.filter.FilteringIterator;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public class SingleLocusContextIterator extends LocusContextIterator {

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------
    private final PushbackIterator<SAMRecord> it;
    private String contig = null;
    private int position = -1;
    private List<SAMRecord> reads = new ArrayList<SAMRecord>(100);
    private List<Integer> offsets = new ArrayList<Integer>(100);

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public SingleLocusContextIterator(final CloseableIterator<SAMRecord> samIterator) {
        FilteringIterator filterIter = new FilteringIterator(samIterator, new locusStreamFilterFunc());
        this.it = new PushbackIterator<SAMRecord>(filterIter);

    }

    /**
     * Class to filter out un-handle-able reads from the stream.  
     */
    class locusStreamFilterFunc implements SamRecordFilter {
        public boolean filterOut(SAMRecord rec) {
            return rec.getCigar().numCigarElements() > 1;
        }
    }

    public Iterator<LocusContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public boolean hasNext() {
        return it.hasNext();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public LocusContext next() {
        position += 1;

        if ( position != -1 ) {
            cleanReads();
            expandReads();
        }

        if ( reads.isEmpty() ) {
            //  the window is empty, we need to jump to the first pos of the first read in the stream
            SAMRecord read = it.next();
            pushRead(read);
            contig = read.getReferenceName();
            position = read.getAlignmentStart() - 1;
            return next();
        }
        else {
            // at this point, window contains all reads covering the pos, we need to return them
            // and the offsets into each read for this loci
            calcOffsetsOfWindow(position);
            return new LocusContext(new GenomeLoc(contig, position), reads, offsets);
        }
    }

    private void pushRead(SAMRecord read) {
        //System.out.printf("  -> Adding read %s %d-%d flags %s%n", read.getReadName(), read.getAlignmentStart(), read.getAlignmentEnd(), Utils.readFlagsAsString(read));
        reads.add(read);
    }

    class KeepReadPFunc implements Predicate<SAMRecord> {
        public boolean apply(SAMRecord read) {
            return position >= read.getAlignmentStart() &&
                    position < read.getAlignmentEnd() &&
                    read.getReferenceName().equals(contig); // should be index for efficiency
        }
    }
    Predicate KeepReadP = new SingleLocusContextIterator.KeepReadPFunc();

    private void calcOffsetsOfWindow(final int position) {
        offsets.clear();
        for ( SAMRecord read : reads ) {
//            def calcOffset( read ):
//                offset = self.pos - read.start
//                return offset
//
//            offsets = map(calcOffset, self.window)
            final int offset = position - read.getAlignmentStart();
            assert(offset < read.getReadLength() );
            offsets.add(offset);
            //System.out.printf("offsets [%d] %s%n", read.getAlignmentStart(), offsets);
        }
    }

    private void cleanReads() {
        // def keepReadP( read ):
        //     return read.chr == chr and pos >= read.start and pos <= read.end
        // self.window = filter( keepReadP, self.window )
        reads = Utils.filter(KeepReadP, reads);
    }

    private void expandReads() {
//        for read in self.rs:
//            #print 'read', read, pos
//            if read.chr == chr and read.start <= pos and read.end >= pos:
//                self.pushRead(read)
//            else:
//                self.rs.unget( read )
//                #self.rs = chain( [read], self.rs )
//                break
        while ( it.hasNext() ) {
            SAMRecord read = it.next();
            if ( KeepReadP.apply( read ) ) {
                pushRead(read);
            }
            else {
                it.pushback(read);
                break;
            }
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
