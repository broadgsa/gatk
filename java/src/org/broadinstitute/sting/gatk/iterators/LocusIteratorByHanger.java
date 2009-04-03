package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import org.broadinstitute.sting.utils.*;

import java.util.List;
import java.util.Iterator;

import org.broadinstitute.sting.utils.RefHanger;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.apache.log4j.Logger;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public class LocusIteratorByHanger extends LocusIterator {

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(LocusIteratorByHanger.class);

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
    boolean justCleared = false;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public LocusIteratorByHanger(final Iterator<SAMRecord> samIterator) {
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

            logger.debug(String.format("  -> %s:", rhanger.loc));
            for ( int j = 0; j < rhanger.size(); j++ ) {
                SAMRecord read = (SAMRecord)rhanger.get(j);
                int offset = (Integer)ohanger.get(j);
                logger.debug(String.format(" %s(%d)=%s", read.getReadName(), offset, read.getReadString().charAt(offset) ));
            }
            logger.debug(String.format("%n"));

        }        
    }

    public void clear() {
        logger.debug(String.format(("clear() called%n")));
        readHanger.clear();
        offsetHanger.clear();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public LocusContext next() {
//        if ( it.hasNext() && ! readHanger.isEmpty() ) {
//            // todo: this needs to be deleted
//            final SAMRecord read = it.peek();
//            GenomeLoc readLoc = Utils.genomicLocationOf(read);
//            System.out.printf("Comparing %s to %s%n", readLoc, readHanger.getLeftLoc());
//            if ( readLoc.compareTo(readHanger.getLeftLoc()) == -1 ) {
//                clear();
//                return next();
//            }
//        }

        if ( ! currentPositionIsFullyCovered() )
            expandWindow(INCREMENT_SIZE);

        if ( DEBUG ) {
            logger.debug(String.format(("in Next:%n")));
            printState();
        }

        RefHanger.Hanger rhanger = readHanger.popLeft();
        RefHanger.Hanger ohanger = offsetHanger.popLeft();

        return new LocusContext(rhanger.loc, rhanger.data, ohanger.data);
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
                logger.debug(String.format("Processing block %s len=%d%n", block, block.getLength()));
            for ( int i = 0; i < block.getLength(); i++ ) {
                GenomeLoc offset = new GenomeLoc(readLoc.getContig(), block.getReferenceStart() + i);
                readHanger.expandingPut(offset, read);
                offsetHanger.expandingPut(offset, block.getReadStart() + i - 1);
                if ( DEBUG )
                    logger.debug(String.format("  # Added %s%n", offset));
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
        if ( ! it.hasNext() )   // if there are no more reads, we are fully covered
            return true;
        else {
            final SAMRecord read = it.peek();
            GenomeLoc readLoc = GenomeLoc.genomicLocationOf(read);
            final boolean coveredP = currentPositionIsFullyCovered(readLoc);
            //System.out.printf("CoverP = %s => %b%n", readLoc, coveredP);
            return coveredP;
        }
    }

    private final void expandWindow(final int incrementSize) {
        if ( DEBUG ) {
            logger.debug(String.format("entering expandWindow..., hasNext=%b%n", it.hasNext()));
            printState();
        }

        while ( it.hasNext() ) {
            if ( DEBUG ) {
                logger.debug(String.format("Expanding window%n"));
                printState();
            }
            
            SAMRecord read = it.next();
            justCleared = false;

            GenomeLoc readLoc = GenomeLoc.genomicLocationOf(read);
            if ( DEBUG ) {
                logger.debug(String.format("  Expanding window sizes %d with %d : left=%s, right=%s, readLoc = %s, cmp=%d%n",
                        readHanger.size(), incrementSize,
                        readHanger.hasHangers() ? readHanger.getLeftLoc() : "NA",
                        readHanger.hasHangers() ? readHanger.getRightLoc() : "NA",
                        readLoc,
                        readHanger.hasHangers() ? readLoc.compareTo(readHanger.getLeftLoc()) : -100));
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

//    private final void expandWindow(final int incrementSize) {
//        if ( DEBUG ) {
//            System.out.printf("entering expandWindow..., hasNext=%b%n", it.hasNext());
//            printState();
//        }
//
//        while ( it.hasNext() ) {
//            if ( DEBUG ) {
//                System.out.printf("Expanding window%n");
//                printState();
//            }
//
//            try {
//                SAMRecord read = it.next();
//                justCleared = false;
//
//                GenomeLoc readLoc = Utils.genomicLocationOf(read);
//                if ( DEBUG ) {
//                    System.out.printf("  Expanding window sizes %d with %d : left=%s, right=%s, readLoc = %s, cmp=%d%n",
//                            readHanger.size(), incrementSize,
//                            readHanger.hasHangers() ? readHanger.getLeftLoc() : "NA",
//                            readHanger.hasHangers() ? readHanger.getRightLoc() : "NA",
//                            readLoc,
//                            readHanger.hasHangers() ? readLoc.compareTo(readHanger.getLeftLoc()) : -100);
//                }
//                //if ( readHanger.size() >= incrementSize ) {
//                //if ( readHanger.hasHangers() && readLoc.compareTo(readHanger.getLeftLoc()) == 1) {
//                if ( readHanger.hasHangers() && readLoc.distance(readHanger.getLeftLoc()) >= incrementSize ) {
//                    // We've collected up enough reads
//                    it.pushback(read);
//                    break;
//                }
//                else
//                    hangRead(read);
//            }
//            catch ( RuntimeIOException rio ) {
//                System.out.printf("Clearing state...");
//                // todo: good god, this is dangerous, we are reseting state because the reads are out of order
//                if ( ! justCleared ) {
//                    rio.printStackTrace();
//                    justCleared = true;
//                    clear();
//                }
//                it.next(); // throw away the offending read
//                expandWindow(INCREMENT_SIZE);
//                return;
//            }
//        }
//    }


    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}