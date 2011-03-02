package org.broadinstitute.sting.utils.sam;

import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
//import org.broadinstitute.sting.utils.SimpleTimer;

import java.io.File;
import java.util.*;

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

/**
 * A locally resorting, mate fixing sam file writer that supports an idiom where reads are only moved around if
 * the ISIZE of the pair is < X and reads are not allowed to move any more than Y bp from their original positions.
 *
 * To understand this data structure, let's begin by asking -- when are we certain we know the position of read R added
 * to the writer and its mate M given that R has been added to the writer (but M may not be), their ISIZE in R, at the
 * moment that a read K is added to the writer, under the constraints X and Y?  Complex I know.  First, because
 * reads cannot move more than Y bp in either direction, we know that R originated at most R.pos + Y bp from its
 * current position.  Also, we know that K is at most K.pos + Y bp from it's original position.  If R is maximally
 * shifted to the right, and K shifted to the left, then they could at most move 2Y together.  So if the distance
 * between R and K > 2Y, we know that there are no reads left in the original stream that could be moved before R.
 *
 * Now, we also need to be certain if we have a mate pair M, that won't emit R before we can incorporate any move of
 * M into the mate pair info R.  There are two cases to consider here:
 *
 * If ISIZE > X, we know that we won't move M when we see it, so we can safely emit R knowing that
 * M is fixed in place.
 *
 * If ISIZE <= X, M might be moved, and it we have to wait until we see M in the stream to know it's position.
 * So R must be buffered until either M arrives, or we see a read K that's more than 2Y units past the original position
 * of M.
 *
 * So the worst-case memory consumption here is proportional to the number of reads
 * occurring between R and M + 2 Y, and so is proportional to the depth of the data and X and Y.
 *
 * This leads to the following simple algorithm:
 *
 * addAlignment(newRead):
 *   addReadToListOfReads(newRead)
 *   update mate pair of newRead if present in list of reads
 *
 *   for ( read in list of reads [in order of increasing read.pos] ):
 *     if read.pos < newRead.pos - 2Y && (read.isize >= X || read.matePos < newRead.pos - 2 * Y):
 *        emit read and remove from list of reads
 *     else:
 *        break
 *
 * @author depristo
 * @version 0.1
 */
public class ConstrainedMateFixingSAMFileWriter implements SAMFileWriter {
    protected static final Logger logger = Logger.getLogger(ConstrainedMateFixingSAMFileWriter.class);
    private static final boolean DEBUG = false;

    /** How often do we check whether we want to emit reads? */
    private final static int EMIT_FREQUENCY = 1000;

    /**
     * How much could a single read move in position from its original position?
     */
    private int MAX_POS_MOVE_ALLOWED;

    /** how we order our SAM records */
    private final SAMRecordComparator comparer = new SAMRecordCoordinateComparator();


    /** The place where we ultimately write out our records */
    final SAMFileWriter finalDestination;

    /**
     * what is the maximum isize of a pair of reads that can move?  Reads with isize > this value
     * are assumes to not be allowed to move in the incoming read stream.
     */
    final int maxInsertSizeForMovingReadPairs;

    int counter = 0;
    int maxReadsInQueue = 0;

    /** read.name -> records */
    HashMap<String, SAMRecord> forMateMatching = new HashMap<String, SAMRecord>();
    Queue<SAMRecord> waitingReads = new PriorityQueue<SAMRecord>(1000, comparer);

    //private SimpleTimer timer = new SimpleTimer("ConstrainedWriter");
    //private long PROGRESS_PRINT_FREQUENCY = 10 * 1000;             // in milliseconds
    //private long lastProgressPrintTime = -1;                       // When was the last time we printed progress log?


    /**
     *
     * @param header
     * @param outputFile
     * @param compressionLevel
     * @param maxInsertSizeForMovingReadPairs
     */
    public ConstrainedMateFixingSAMFileWriter(final SAMFileHeader header,
                                              final File outputFile,
                                              final int compressionLevel,
                                              final int maxInsertSizeForMovingReadPairs,
                                              final int maxMoveAllowed) {
        this(new SAMFileWriterFactory().makeBAMWriter(header, true, outputFile, compressionLevel),
                maxInsertSizeForMovingReadPairs,
                maxMoveAllowed);
    }

    public ConstrainedMateFixingSAMFileWriter(final SAMFileWriter finalDestination,
                                              final int maxInsertSizeForMovingReadPairs,
                                              final int maxmoveAllowed) {
        this.finalDestination = finalDestination;
        this.maxInsertSizeForMovingReadPairs = maxInsertSizeForMovingReadPairs;
        this.MAX_POS_MOVE_ALLOWED = maxmoveAllowed;

        //timer.start();
        //lastProgressPrintTime = timer.currentTime();
    }

    public int getMaxReadsInQueue() { return maxReadsInQueue; }
    public int getNReadsInQueue() { return waitingReads.size(); }

    /**
     * Retrieves the header to use when creating the new SAM file.
     * @return header to use when creating the new SAM file.
     */
    public SAMFileHeader getFileHeader() {
        return finalDestination.getFileHeader();
    }

    private boolean noReadCanMoveBefore(int pos, SAMRecord addedRead) {
        return pos + 2 * MAX_POS_MOVE_ALLOWED < addedRead.getAlignmentStart();
    }


    /**
     * @{inheritDoc}
     */
    public void addAlignment( SAMRecord newRead ) {
        if ( DEBUG ) logger.info("New read pos " + newRead.getAlignmentStart() + " OP = " + newRead.getAttribute("OP"));

        //final long curTime = timer.currentTime();
        //if ( curTime - lastProgressPrintTime > PROGRESS_PRINT_FREQUENCY ) {
        //    lastProgressPrintTime = curTime;
        //    System.out.println("WaitingReads.size = " + waitingReads.size() + ", forMateMatching.size = " + forMateMatching.size());
        //}

        // if the new read is on a different contig, then we need to flush the queue and clear the map
        if ( waitingReads.size() > 0 && waitingReads.peek().getReferenceIndex() != newRead.getReferenceIndex()) {
            if ( DEBUG ) logger.warn("Flushing queue on move to new contig: " + newRead.getReferenceName());

            while ( ! waitingReads.isEmpty() ) {
                // emit to disk
                finalDestination.addAlignment(waitingReads.remove());
            }

            forMateMatching.clear();
        }

        // fix mates, as needed
        // Since setMateInfo can move reads, we potentially need to remove the mate, and requeue
        // it to ensure proper sorting
        if ( newRead.getReadPairedFlag() ) {
            SAMRecord mate = forMateMatching.get(newRead.getReadName());
            if ( mate != null ) {
                // Frustratingly, Picard's setMateInfo() method unaligns (by setting the reference contig
                // to '*') read pairs when both of their flags have the unmapped bit set.  This is problematic
                // when trying to emit reads in coordinate order because all of a sudden we have reads in the
                // middle of the bam file that now belong at the end - and any mapped reads that get emitted
                // after them trigger an exception in the writer.  For our purposes, because we shouldn't be
                // moving read pairs when they are both unmapped anyways, we'll just not run fix mates on them.
                boolean doNotFixMates = newRead.getReadUnmappedFlag() && mate.getReadUnmappedFlag();
                if ( !doNotFixMates ) {

                    boolean reQueueMate = mate.getReadUnmappedFlag() && ! newRead.getReadUnmappedFlag();
                    if ( reQueueMate ) {
                        // the mate was unmapped, but newRead was mapped, so the mate may have been moved
                        // to be next-to newRead, so needs to be reinserted into the waitingReads queue
                        // note -- this must be called before the setMateInfo call below
                        if ( ! waitingReads.remove(mate) )
                            throw new ReviewedStingException("BUG: removal of mate failed at " + mate);
                    }

                    // we've already seen our mate -- set the mate info and remove it from the map
                    SamPairUtil.setMateInfo(mate, newRead, null);
                    if ( reQueueMate ) waitingReads.add(mate);
                }

                forMateMatching.remove(newRead.getReadName());
            } else {
                forMateMatching.put(newRead.getReadName(), newRead);
            }
        }

        waitingReads.add(newRead);
        maxReadsInQueue = Math.max(maxReadsInQueue, waitingReads.size());

        if ( ++counter % EMIT_FREQUENCY == 0 ) {
            while ( ! waitingReads.isEmpty() ) { // there's something in the queue
                SAMRecord read = waitingReads.peek();

                if ( noReadCanMoveBefore(read.getAlignmentStart(), newRead) &&
                        (iSizeTooBigToMove(read)                                           // we won't try to move such a read
                                || ! read.getReadPairedFlag()                                     // we're not a paired read
                                || read.getReadUnmappedFlag() && read.getMateUnmappedFlag()       // both reads are unmapped
                                || noReadCanMoveBefore(read.getMateAlignmentStart(), newRead ) ) ) { // we're already past where the mate started

                    // remove reads from the map that we have emitted -- useful for case where the mate never showed up
                    forMateMatching.remove(read.getReadName());

                    if ( DEBUG )
                        logger.warn(String.format("EMIT!  At %d: read %s at %d with isize %d, mate start %d, op = %s",
                                newRead.getAlignmentStart(), read.getReadName(), read.getAlignmentStart(),
                                read.getInferredInsertSize(), read.getMateAlignmentStart(), read.getAttribute("OP")));
                    // emit to disk
                    finalDestination.addAlignment(waitingReads.remove());
                } else {
                    if ( DEBUG )
                        logger.warn(String.format("At %d: read %s at %d with isize %d couldn't be emited, mate start %d",
                                newRead.getAlignmentStart(), read.getReadName(), read.getAlignmentStart(), read.getInferredInsertSize(), read.getMateAlignmentStart()));
                    break;
                }
            }

            if ( DEBUG ) logger.warn(String.format("At %d: Done with emit cycle", newRead.getAlignmentStart()));
        }
    }

    /**
     * Returns true if the read shouldn't be moved given the constraints of this SAMFileWriter
     * @param read
     * @return
     */
    public boolean iSizeTooBigToMove(SAMRecord read) {
        return iSizeTooBigToMove(read, maxInsertSizeForMovingReadPairs);               // we won't try to move such a read
    }

    public static boolean iSizeTooBigToMove(SAMRecord read, int maxInsertSizeForMovingReadPairs) {
        return ( read.getReadPairedFlag() && ! read.getMateUnmappedFlag() && read.getReferenceName() != read.getMateReferenceName() ) // maps to different chromosomes
                || Math.abs(read.getInferredInsertSize()) > maxInsertSizeForMovingReadPairs;     // we won't try to move such a read
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() ) { // there's something in the queue
            finalDestination.addAlignment(waitingReads.remove());
        }
        finalDestination.close();
    }
}
