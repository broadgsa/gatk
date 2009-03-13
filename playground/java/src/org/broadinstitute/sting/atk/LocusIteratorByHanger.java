package org.broadinstitute.sting.atk;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import org.broadinstitute.sting.utils.*;

import java.util.List;
import java.util.Iterator;

import org.broadinstitute.sting.utils.RefHanger;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public class LocusIteratorByHanger extends LocusIterator {

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------
    private final PushbackIterator<SAMRecord> it;

    private RefHanger<SAMRecord> readHanger = new RefHanger<SAMRecord>();
    private RefHanger<Integer> offsetHanger = new RefHanger<Integer>();
    final int INCREMENT_SIZE = 100;
    final boolean DEBUG = false;

    /**
     * Useful class for forwarding on locusContext data from this iterator
     */
    public class MyLocusContext implements LocusContext {
        GenomeLoc loc = null;
        private List<SAMRecord> reads = null;
        private List<Integer> offsets = null;

        private MyLocusContext(GenomeLoc loc, List<SAMRecord> reads, List<Integer> offsets) {
            this.loc = loc;
            this.reads = reads;
            this.offsets = offsets;
        }

        public String getContig() { return getLocation().getContig(); }
        public long getPosition() { return getLocation().getStart(); }
        public GenomeLoc getLocation() { return loc; }

        public List<SAMRecord> getReads() { return reads; }
        public List<Integer> getOffsets() { return offsets; }
        public int numReads() { return reads.size(); }
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public LocusIteratorByHanger(final CloseableIterator<SAMRecord> samIterator) {
        this.it = new PushbackIterator<SAMRecord>(samIterator);
    }

    public Iterator<LocusContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public boolean hasNext() {
        return readHanger.hasHangers() || it.hasNext();
    }

    public void printState() {
        for ( int i = 0; i < readHanger.size(); i++ ) {
            RefHanger.Hanger rhanger = readHanger.getHanger(i);
            RefHanger.Hanger ohanger = offsetHanger.getHanger(i);

            System.out.printf("  -> %s:", rhanger.loc);
            for ( int j = 0; j < rhanger.size(); j++ ) {
                SAMRecord read = (SAMRecord)rhanger.get(j);
                int offset = (Integer)ohanger.get(j);
                System.out.printf(" %s(%d)=%s", read.getReadName(), offset, read.getReadString().charAt(offset) );
            }
            System.out.printf("%n");

        }        
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public MyLocusContext next() {
        if ( ! currentPositionIsFullyCovered() )
            expandWindow(INCREMENT_SIZE);

        if ( DEBUG ) {
            System.out.printf("in Next:%n");
            printState();
        }

        RefHanger.Hanger rhanger = readHanger.popLeft();
        RefHanger.Hanger ohanger = offsetHanger.popLeft();

        return new MyLocusContext(rhanger.loc, rhanger.data, ohanger.data);
    }

    protected void hangRead(final SAMRecord read) {
        GenomeLoc readLoc = new GenomeLoc(read.getReferenceName(), read.getAlignmentStart());
        //System.out.printf("Adding read %s at %d%n", read.getReadName(), read.getAlignmentStart());
        /*
        for ( int i = 0; i < read.getReadLength(); i++ ) {
            GenomeLoc offset = new GenomeLoc(readLoc.getContig(), readLoc.getStart() + i);
            readHanger.expandingPut(offset, read);
            offsetHanger.expandingPut(offset, i);
        }
        */

        for ( AlignmentBlock block : read.getAlignmentBlocks() ) {
            if ( DEBUG )
                System.out.printf("Processing block %s len=%d%n", block, block.getLength());
            for ( int i = 0; i < block.getLength(); i++ ) {
                GenomeLoc offset = new GenomeLoc(readLoc.getContig(), block.getReferenceStart() + i);
                readHanger.expandingPut(offset, read);
                offsetHanger.expandingPut(offset, block.getReadStart() + i - 1);
                if ( DEBUG )
                    System.out.printf("  # Added %s%n", offset);
            }
        }
    }

    private final boolean currentPositionIsFullyCovered(final GenomeLoc nextReadInStreamLoc) {
        if ( readHanger.isEmpty() )
            // If the buffer is empty, we're definitely not done
            return false;

        if ( nextReadInStreamLoc.compareTo(readHanger.getLeftLoc()) == 1 )
            // the next read in the stream is beyond the left most position, so it's fully covered
            return true;
        else
            // not fully covered
            return false;
    }

    private final boolean currentPositionIsFullyCovered() {
        final SAMRecord read = it.next();
        GenomeLoc readLoc = new GenomeLoc(read.getReferenceName(), read.getAlignmentStart());
        final boolean coveredP = currentPositionIsFullyCovered(readLoc);
        if ( coveredP )
            it.pushback(read);
        return coveredP;
    }

    private final void expandWindow(final int incrementSize) {
        while ( it.hasNext() ) {
            if ( DEBUG ) {
                System.out.printf("Expanding window%n");
                printState();
            }
            
            SAMRecord read = it.next();

            GenomeLoc readLoc = new GenomeLoc(read.getReferenceName(), read.getAlignmentStart());
            if ( DEBUG ) {
                System.out.printf("  Expanding window sizes %d with %d : left=%s, right=%s, readLoc = %s, cmp=%d%n",
                    readHanger.size(), incrementSize,
                    readHanger.hasHangers() ? readHanger.getLeftLoc() : "NA", 
                    readHanger.hasHangers() ? readHanger.getRightLoc() : "NA",
                    readLoc,
                    readHanger.hasHangers() ? readLoc.compareTo(readHanger.getLeftLoc()) : -100);
            }
            //if ( readHanger.size() >= incrementSize ) {
            //if ( readHanger.hasHangers() && readLoc.compareTo(readHanger.getLeftLoc()) == 1) {
            if ( readHanger.hasHangers() && readLoc.distance(readHanger.getLeftLoc()) >= incrementSize ) {
                // We've collected up enough reads
                it.pushback(read);
                break;
            }
            else
                hangRead(read);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}