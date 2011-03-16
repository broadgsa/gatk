/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.Bin;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKChunk;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Writes schedules for a single BAM file to a target output file.
 */
public class BAMSchedule implements CloseableIterator<BAMScheduleEntry> {
    /**
     * File in which to store schedule data.
     */
    private final File scheduleFile;

    /**
     * File channel for the schedule file.
     */
    private final FileChannel scheduleFileChannel;

    /**
     * Next schedule entry to be returned.  Null if no additional entries are present.
     */
    private BAMScheduleEntry nextScheduleEntry;

    /**
     * Sizes of ints and longs in bytes.
     */
    private static final int INT_SIZE_IN_BYTES = Integer.SIZE / 8;
    private static final int LONG_SIZE_IN_BYTES = Long.SIZE / 8;    

    /**
     * Create a new BAM schedule based on the given index.
     * @param index index to convert to a schedule.
     */
    public BAMSchedule(final GATKBAMIndex index, final int referenceSequence) {
        try {
            scheduleFile = File.createTempFile("bamschedule."+referenceSequence,null);
            scheduleFileChannel = new RandomAccessFile(scheduleFile,"rw").getChannel();
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to create BAM schedule file.",ex);
        }
        scheduleFile.deleteOnExit();

        int currentBinInLowestLevel = GATKBAMIndex.getFirstBinInLevel(GATKBAMIndex.getNumIndexLevels()-1) - 1;
        while(++currentBinInLowestLevel < GATKBAMIndex.MAX_BINS) {
            BAMScheduleEntry scheduleEntry = BAMScheduleEntry.query(index,referenceSequence,currentBinInLowestLevel);
            if(scheduleEntry.fileSpan.isEmpty())
                continue;

            // File format is binary in little endian; start of region, end of region, num chunks, then the chunks themselves.
            ByteBuffer buffer = ByteBuffer.allocateDirect(2*INT_SIZE_IN_BYTES + INT_SIZE_IN_BYTES + scheduleEntry.fileSpan.getGATKChunks().size()*LONG_SIZE_IN_BYTES*2);
            buffer.order(ByteOrder.LITTLE_ENDIAN);
            buffer.putInt(scheduleEntry.start);
            buffer.putInt(scheduleEntry.stop);
            buffer.putInt(scheduleEntry.fileSpan.getGATKChunks().size());
            for(GATKChunk chunk: scheduleEntry.fileSpan.getGATKChunks()) {
                buffer.putLong(chunk.getChunkStart());
                buffer.putLong(chunk.getChunkEnd());
            }

            // Prepare buffer for writing
            buffer.flip();

            try {
                scheduleFileChannel.write(buffer);
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to create BAM schedule file.",ex);
            }
        }

        // Move file pointer back to the start.
        try {
            scheduleFileChannel.position(0L);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to rewind BAM schedule file.",ex);
        }

        advance();
    }

    /**
     * Determine whether more ScheduleEntries are present in the iterator.
     * @return Next schedule entry to parse.
     */
    @Override
    public boolean hasNext() {
        return nextScheduleEntry != null;    
    }

    /**
     * Retrieve the next schedule entry in the list.
     * @return next schedule entry in the queue.
     */
    @Override
    public BAMScheduleEntry next() {
        BAMScheduleEntry currentScheduleEntry = nextScheduleEntry;
        advance();
        return currentScheduleEntry;
    }

    /**
     * Close down and delete the file.
     */
    @Override
    public void close() {
        try {
            scheduleFileChannel.close();
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to close schedule file.");
        }
    }

    /**
     * Advance to the next schedule entry.
     */
    private void advance() {
        nextScheduleEntry = null;

        ByteBuffer buffer = ByteBuffer.allocateDirect(2*INT_SIZE_IN_BYTES+INT_SIZE_IN_BYTES);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        int results;

        try {
            results = scheduleFileChannel.read(buffer);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to read start, stop and chunk sizes from schedule file channel.",ex);
        }

        // No more data to read.
        if(results <= 0)
            return;

        // Reorient buffer for reading.
        buffer.flip();

        final int start = buffer.getInt();
        final int stop = buffer.getInt();
        final int numChunks = buffer.getInt();

        GATKChunk[] chunks = new GATKChunk[numChunks];
        buffer = ByteBuffer.allocateDirect(numChunks * 2 * LONG_SIZE_IN_BYTES);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        
        try {
            scheduleFileChannel.read(buffer);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to read chunk data from schedule file channel.",ex);
        }

        // Reposition for reading.
        buffer.flip();

        // Read out chunk data.
        for(int i = 0; i < numChunks; i++)
            chunks[i] = new GATKChunk(buffer.getLong(),buffer.getLong());

        // Prep the iterator for the next schedule entry.
        nextScheduleEntry = new BAMScheduleEntry(start,stop,new GATKBAMFileSpan(chunks));
    }

    @Override
    public void remove() { throw new UnsupportedOperationException("Unable to remove from a schedule iterator."); }
}

/**
 * A single proto-shard to be processed.
 */
class BAMScheduleEntry {
    /**
     * Starting position for the genomic entry.
     */
    public final int start;

    /**
     * Ending position for the genomic entry.
     */
    public final int stop;

    /**
     * The spans representing the given region.
     */
    public final GATKBAMFileSpan fileSpan;

    BAMScheduleEntry(final int start, final int stop, final GATKBAMFileSpan fileSpan) {
        this.start = start;
        this.stop = stop;
        this.fileSpan = fileSpan;
    }

    public static BAMScheduleEntry query(final GATKBAMIndex index, final int referenceSequence, final int binNumber) {
        final Bin bin = new Bin(referenceSequence,binNumber);
        final int start = index.getFirstLocusInBin(bin);
        final int stop = index.getLastLocusInBin(bin);
        final GATKBAMFileSpan fileSpan = index.getSpanOverlapping(bin);
        return new BAMScheduleEntry(start,stop,fileSpan);
    }

    /**
     * Returns true if the location of this bin tree is before the given position.
     * @param locus Locus to test.
     * @return True if this bin sits completely before the given locus; false otherwise.
     */
    public boolean isBefore(final GenomeLoc locus) {
        return stop < locus.getStart();
    }

    /**
     * Checks overlap between this bin tree and other bin trees.
     * @param position the position over which to detect overlap.
     * @return True if the segment overlaps.  False otherwise.
     */
    public boolean overlaps(final GenomeLoc position) {
        return !(position.getStop() < start || position.getStart() > stop);
    }
}
