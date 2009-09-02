/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.RefHanger;
import org.broadinstitute.sting.utils.Utils;

import java.util.Iterator;

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
    private Reads readInfo;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public LocusIteratorByHanger(final Iterator<SAMRecord> samIterator, Reads readInformation) {
        this.it = new PushbackIterator<SAMRecord>(samIterator);
        this.readInfo = readInformation;
    }

    public Iterator<AlignmentContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public boolean hasNext() {
        return readHanger.hasHangers() || it.hasNext();
    }

    public void printState() {
        for (int i = 0; i < readHanger.size(); i++) {
            RefHanger.Hanger rhanger = readHanger.getHanger(i);
            RefHanger.Hanger ohanger = offsetHanger.getHanger(i);

            logger.debug(String.format("printState(): location %s", rhanger.loc));
            for (int j = 0; j < rhanger.size(); j++) {
                SAMRecord read = (SAMRecord) rhanger.get(j);
                int offset = (Integer) ohanger.get(j);
                logger.debug(String.format("  read: %s(%d)=%s", read.getReadName(), offset, read.getReadString().charAt(offset)));
            }
        }
    }

    public void clear() {
        logger.debug(String.format(("clear() called")));
        readHanger.clear();
        offsetHanger.clear();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public AlignmentContext next() {
        if (!currentPositionIsFullyCovered())
            expandWindow(INCREMENT_SIZE, readInfo.getMaxReadsAtLocus());

        if (DEBUG) {
            logger.debug("in Next:");
            printState();
        }

        RefHanger.Hanger rhanger = readHanger.popLeft();
        RefHanger.Hanger ohanger = offsetHanger.popLeft();

        if (rhanger.size() == 0) {
            // we are in the case where there are no reads (we're inside a large indel without any reads
            // so recursively call next.  This is safe because having a position without any reads
            // implies that there are positions to the right with reads
            //System.out.printf("***** skipping reads%n");
            return next();
        } else {
            return new AlignmentContext(rhanger.loc, rhanger.data, ohanger.data);
        }
    }

    protected boolean hangRead(final SAMRecord read, final int maximumPileupSize, boolean warned) {
        GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);

        for (AlignmentBlock block : read.getAlignmentBlocks()) {
            if (DEBUG) logger.debug(String.format("Processing block %s len=%d", block, block.getLength()));
            for (int i = 0; i < block.getLength(); i++) {
                // check to see if we've exceeded the maximum number of reads in the pile-up
                GenomeLoc offset = GenomeLocParser.createGenomeLoc(readLoc.getContigIndex(), block.getReferenceStart() + i);
                int hangerSize = (readHanger.hasLocation(offset)) ? readHanger.getHanger(offset).size() : -1;
                if (hangerSize >= maximumPileupSize) {
                    if (!warned) {
                        warned = true;
                        Utils.warnUser("Unable to add a read, we're over the hanger limit of " + maximumPileupSize + " at location " + readHanger.getLeftLoc());
                    }
                } else {
                    readHanger.expandingPut(offset, read);
                    offsetHanger.expandingPut(offset, block.getReadStart() + i - 1);
                }                    
                if (DEBUG) logger.debug(String.format("  # Added %s", offset));
            }
        }
        return warned; // did we warn the user about this location?
    }

    private final boolean currentPositionIsFullyCovered(final GenomeLoc nextReadInStreamLoc) {
        if (readHanger.isEmpty())
            // If the buffer is empty, we're definitely not done
            return false;

        if (nextReadInStreamLoc.compareTo(readHanger.getLeftLoc()) == 1)
            // the next read in the stream is beyond the left most position, so it's fully covered
            return true;
        else
            // not fully covered
            return false;
    }

    private final boolean currentPositionIsFullyCovered() {
        if (!it.hasNext())   // if there are no more reads, we are fully covered
            return true;
        else {
            final SAMRecord read = it.peek();
            GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
            final boolean coveredP = currentPositionIsFullyCovered(readLoc);
            //System.out.printf("CoverP = %s => %b%n", readLoc, coveredP);
            return coveredP;
        }
    }

    private final void expandWindow(final int incrementSize, final int maximumPileupSize) {
        // todo: reenable
        if (false && incrementSize != 1) {
            Utils.scareUser(String.format("BUG: IncrementSize=%d != 1, the codebase doesn't support this extension strategy yet", incrementSize));
        }

        if (DEBUG) {
            logger.debug(String.format("entering expandWindow..., hasNext=%b", it.hasNext()));
            printState();
        }
        Boolean warned = false; // warn them once per locus
        while (it.hasNext()) {
            if (DEBUG) {
                logger.debug(String.format("Expanding window"));
                printState();
            }

            SAMRecord read = it.next();
            justCleared = false;

            GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
            if (DEBUG) {
                logger.debug(String.format("  Expanding window sizes %d with %d : left=%s, right=%s, readLoc = %s, cmp=%d",
                        readHanger.size(), incrementSize,
                        readHanger.hasHangers() ? readHanger.getLeftLoc() : "NA",
                        readHanger.hasHangers() ? readHanger.getRightLoc() : "NA",
                        readLoc,
                        readHanger.hasHangers() ? readLoc.compareTo(readHanger.getLeftLoc()) : -100));
            }
            //if ( readHanger.size() >= incrementSize ) {
            //if ( readHanger.hasHangers() && readLoc.compareTo(readHanger.getLeftLoc()) == 1) {
            if (readHanger.hasHangers() && readLoc.distance(readHanger.getLeftLoc()) >= incrementSize) {
                // We've collected up enough reads
                it.pushback(read);
                break;
            } else
                warned = hangRead(read,maximumPileupSize,warned);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
