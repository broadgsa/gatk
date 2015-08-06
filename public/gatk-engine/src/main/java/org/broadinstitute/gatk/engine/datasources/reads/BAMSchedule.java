/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.Bin;
import htsjdk.samtools.GATKBAMFileSpan;
import htsjdk.samtools.GATKChunk;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * Writes schedules for a single BAM file to a target output file.
 */
public class BAMSchedule implements CloseableIterator<BAMScheduleEntry> {
    /**
     * File in which to store schedule data.
     */
    private File scheduleFile;

    /**
     * File channel for the schedule file.
     */
    private FileChannel scheduleFileChannel;

    /**
     * The definitive, sorted list of reader IDs.  Order is important here: the order
     * in which the reader IDs are presented here maps to the order in which they appear in the file. 
     */
    private final List<SAMReaderID> readerIDs = new ArrayList<SAMReaderID>();

    /**
     * Iterators over the schedule.  Stored in the same order as readerIDs, above.
     */
    private final List<PeekableIterator<BAMScheduleEntry>> scheduleIterators = new ArrayList<PeekableIterator<BAMScheduleEntry>>();

    /**
     * Next schedule entry to be returned.  Null if no additional entries are present.
     */
    private BAMScheduleEntry nextScheduleEntry;

    /**
     * Reference sequence for which to write the schedule.
     */
    private final int referenceSequence;

    /**
     * Sizes of ints and longs in bytes.
     */
    private static final int INT_SIZE_IN_BYTES = Integer.SIZE / 8;
    private static final int LONG_SIZE_IN_BYTES = Long.SIZE / 8;    

    /**
     * Create a new BAM schedule based on the given index.
     * @param dataSource The SAM data source to use.
     * @param intervals List of 
     */
    public BAMSchedule(final SAMDataSource dataSource, final List<GenomeLoc> intervals) {
        if(intervals.isEmpty())
            throw new ReviewedGATKException("Tried to write schedule for empty interval list.");

        referenceSequence = dataSource.getHeader().getSequence(intervals.get(0).getContig()).getSequenceIndex();

        createScheduleFile();

        readerIDs.addAll(dataSource.getReaderIDs());

        for(final SAMReaderID reader: readerIDs) {
            final GATKBAMIndex index = dataSource.getIndex(reader);
            final GATKBAMIndexData indexData = index.readReferenceSequence(referenceSequence);

            int currentBinInLowestLevel = GATKBAMIndex.getFirstBinInLevel(GATKBAMIndex.getNumIndexLevels()-1);
            Iterator<GenomeLoc> locusIterator = intervals.iterator();
            GenomeLoc currentLocus = locusIterator.next();

            final long readerStartOffset = position();

            int maxChunkCount = 0;

            while(currentBinInLowestLevel < GATKBAMIndex.MAX_BINS && currentLocus != null) {
                final Bin bin = new Bin(referenceSequence,currentBinInLowestLevel);
                final int binStart = index.getFirstLocusInBin(bin);
                final int binStop = index.getLastLocusInBin(bin);

                // In required, pull bin iterator ahead to the point of the next GenomeLoc.
                if(binStop < currentLocus.getStart()) {
                    currentBinInLowestLevel++;
                    continue;
                }

                // At this point, the bin stop is guaranteed to be >= the start of the locus.
                // If the bins have gone past the current locus, update the current locus if at all possible.
                if(binStart > currentLocus.getStop()) {
                    currentLocus = locusIterator.hasNext() ? locusIterator.next() : null;
                    continue;
                }

                // Code at this point knows that the current bin is neither before nor after the current locus,
                // so it must overlap.  Add this region to the filesystem.
                final GATKBAMFileSpan fileSpan = indexData.getSpanOverlapping(bin);

                if(!fileSpan.isEmpty()) {
                    // File format is binary in little endian; start of region, end of region, num chunks, then the chunks themselves.
                    ByteBuffer buffer = allocateByteBuffer(2*INT_SIZE_IN_BYTES + INT_SIZE_IN_BYTES + fileSpan.getGATKChunks().size()*LONG_SIZE_IN_BYTES*2);
                    buffer.putInt(binStart);
                    buffer.putInt(binStop);
                    buffer.putInt(fileSpan.getGATKChunks().size());
                    for(GATKChunk chunk: fileSpan.getGATKChunks()) {
                        buffer.putLong(chunk.getChunkStart());
                        buffer.putLong(chunk.getChunkEnd());
                    }
                    maxChunkCount = Math.max(maxChunkCount,fileSpan.getGATKChunks().size());

                    // Prepare buffer for writing
                    buffer.flip();

                    // And write.
                    write(buffer);
                }

                currentBinInLowestLevel++;
            }

            final long readerStopOffset = position();

            scheduleIterators.add(new PeekableIterator<BAMScheduleEntry>(new BAMScheduleIterator(reader,readerStartOffset,readerStopOffset,maxChunkCount)));

            // Iterator initialization might move the file pointer.  Make sure it gets reset back to where it was before iterator initialization.
            position(readerStopOffset);
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
            throw makeIOFailureException(true, "Unable to close schedule file.", ex);
        }
    }

    /**
     * Convenience routine for creating UserExceptions
     * @param wasWriting
     * @param message
     * @param e
     * @return
     */
    private final GATKException makeIOFailureException(final boolean wasWriting, final String message, final Exception e) {
        if ( wasWriting ) {
            if ( e == null )
                return new UserException.CouldNotCreateOutputFile(scheduleFile, message);
            else
                return new UserException.CouldNotCreateOutputFile(scheduleFile, message, e);
        } else {
            if ( e == null )
                return new UserException.CouldNotReadInputFile(scheduleFile, message);
            else
                return new UserException.CouldNotReadInputFile(scheduleFile, message, e);
        }
    }

    /**
     * Advance to the next schedule entry.
     */
    private void advance() {
        nextScheduleEntry = null;

        BitSet selectedIterators = new BitSet(readerIDs.size());
        int currentStart = Integer.MAX_VALUE;
        int currentStop = Integer.MAX_VALUE;

        // Select every iterator whose next element is the lowest element in the list.
        for(int reader = 0; reader < scheduleIterators.size(); reader++) {
            PeekableIterator<BAMScheduleEntry> scheduleIterator = scheduleIterators.get(reader);
            if(!scheduleIterator.hasNext())
                continue;

            // If the iterator starts after this one, skip over it.
            if(scheduleIterator.peek().start > currentStart)
                continue;

            // If the iterator starts at the same point as this one, add it to the list.
            if(scheduleIterator.peek().start == currentStart) {
                selectedIterators.set(reader);
                currentStop = Math.min(scheduleIterator.peek().stop,currentStop);
                continue;
            }

            // If the iterator is less than anything seen before it, purge the selections and make this one current.
            if(scheduleIterator.peek().start < currentStart) {
                selectedIterators.clear();
                selectedIterators.set(reader);
                currentStart = scheduleIterator.peek().start;
                currentStop = scheduleIterator.peek().stop;
            }
        }

        // Out of iterators?  Abort early.
        if(selectedIterators.isEmpty())
            return;

        // Create the target schedule entry
        BAMScheduleEntry mergedScheduleEntry = new BAMScheduleEntry(currentStart,currentStop);

        // For each schedule entry with data, load the data into the merged schedule.
        for (int reader = selectedIterators.nextSetBit(0); reader >= 0; reader = selectedIterators.nextSetBit(reader+1)) {
            PeekableIterator<BAMScheduleEntry> scheduleIterator = scheduleIterators.get(reader);
            BAMScheduleEntry individualScheduleEntry = scheduleIterator.peek();
            mergedScheduleEntry.mergeInto(individualScheduleEntry);

            // If the schedule iterator ends after this entry, consume it.
            if(individualScheduleEntry.stop <= currentStop)
                scheduleIterator.next();
        }

        // For each schedule entry without data, add a blank entry.
        for (int reader = selectedIterators.nextClearBit(0); reader < readerIDs.size(); reader = selectedIterators.nextClearBit(reader+1)) {
            mergedScheduleEntry.addFileSpan(readerIDs.get(reader),new GATKBAMFileSpan());
        }

        nextScheduleEntry = mergedScheduleEntry;
    }

    @Override
    public void remove() { throw new UnsupportedOperationException("Unable to remove from a schedule iterator."); }

    /**
     * Create a new schedule file, containing schedule information for all BAM files being dynamically merged.
     */
    private void createScheduleFile() {
        try {
            scheduleFile = File.createTempFile("bamschedule."+referenceSequence,null);
            scheduleFileChannel = new RandomAccessFile(scheduleFile,"rw").getChannel();
        }
        catch(IOException ex) {
            throw new UserException("Unable to create a temporary BAM schedule file.  Please make sure Java can write to the default temp directory or use -Djava.io.tmpdir= to instruct it to use a different temp directory instead.",ex);
        }
        scheduleFile.deleteOnExit();

    }

    /**
     * Creates a new byte buffer of the given size.
     * @param size the size of buffer to allocate.
     * @return Newly allocated byte buffer.
     */
    private ByteBuffer allocateByteBuffer(final int size) {
        ByteBuffer buffer = ByteBuffer.allocate(size);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        return buffer;
    }

    /**
     * Reads the contents at the current position on disk into the given buffer.
     * @param buffer buffer to fill.
     */
    private int read(final ByteBuffer buffer) {
        try {
            return scheduleFileChannel.read(buffer);
        }
        catch(IOException ex) {
            throw makeIOFailureException(false, "Unable to read data from BAM schedule file.", ex);
        }
    }

    private void write(final ByteBuffer buffer) {
        try {
            scheduleFileChannel.write(buffer);
            if(buffer.remaining() > 0)
                throw makeIOFailureException(true, "Unable to write entire buffer to file.", null);
        }
        catch(IOException ex) {
            throw makeIOFailureException(true, "Unable to write data to BAM schedule file.", ex);
        }
    }

    /**
     * Reads the current position from the file channel.
     * @return Current position within file channel.
     */
    private long position() {
        try {
            return scheduleFileChannel.position();
        }
        catch(IOException ex) {
            throw makeIOFailureException(false, "Unable to retrieve position of BAM schedule file.", ex);
        }
    }

    /**
     * Reposition the file channel to the specified offset wrt the start of the file.
     * @param position The position.
     */
    private void position(final long position) {
        try {
            scheduleFileChannel.position(position);
        }
        catch(IOException ex) {
            throw makeIOFailureException(false, "Unable to position BAM schedule file.",ex);
        }
    }

    /**
     * An iterator over the schedule for a single BAM file.
     */
    private class BAMScheduleIterator implements Iterator<BAMScheduleEntry> {
        /**
         * ID of the reader associated with the given schedule.
         */
        private final SAMReaderID reader;

        /**
         * Current position in the file.
         */
        private long currentPosition;

        /**
         * Stopping file position of last bin in file for this reader, exclusive.
         */
        private final long stopPosition;

        /**
         * Byte buffer used to store BAM header info.
         */
        private final ByteBuffer binHeader;

        /**
         * Byte buffer used to store chunk data.
         */
        private final ByteBuffer chunkData;

        public BAMScheduleIterator(final SAMReaderID reader, final long startPosition, final long stopPosition, final int maxChunkCount) {
            this.reader = reader;
            this.currentPosition = startPosition;
            this.stopPosition = stopPosition;
            binHeader = allocateByteBuffer(INT_SIZE_IN_BYTES*3);
            chunkData = allocateByteBuffer(maxChunkCount*LONG_SIZE_IN_BYTES*2);
        }

        @Override
        public boolean hasNext() {
            return currentPosition < stopPosition;
        }

        @Override
        public BAMScheduleEntry next() {
            position(currentPosition);

            // Read data.
            int binHeaderBytesRead = read(binHeader);

            // Make sure we read in a complete bin header:
            if ( binHeaderBytesRead < INT_SIZE_IN_BYTES * 3 ) {
                throw new ReviewedGATKException(String.format("Unable to read a complete bin header from BAM schedule file %s for BAM file %s. " +
                                                               "The BAM schedule file is likely incomplete/corrupt.",
                                                               scheduleFile.getAbsolutePath(), reader.getSamFilePath()));
            }

            // Decode contents.
            binHeader.flip();
            final int start = binHeader.getInt();
            final int stop = binHeader.getInt();
            final int numChunks = binHeader.getInt();

            // Prepare bin buffer for next read.
            binHeader.flip();

            // Prepare a target buffer for chunks.
            GATKChunk[] chunks = new GATKChunk[numChunks];

            // Read all chunk data.
            chunkData.limit(numChunks*LONG_SIZE_IN_BYTES*2);
            long bytesRead = read(chunkData);
            if(bytesRead != numChunks*LONG_SIZE_IN_BYTES*2)
                throw new ReviewedGATKException("Unable to read all chunks from file");

            // Prepare for reading.
            chunkData.flip();

            for(int i = 0; i < numChunks; i++)
                chunks[i] = new GATKChunk(chunkData.getLong(),chunkData.getLong());

            // Prepare chunk buffer for next read.
            chunkData.flip();

            BAMScheduleEntry nextScheduleEntry = new BAMScheduleEntry(start,stop);
            nextScheduleEntry.addFileSpan(reader,new GATKBAMFileSpan(chunks));

            // Reset the position of the iterator at the next contig.
            currentPosition = position();

            return nextScheduleEntry;
        }

        /**
         * Not supported.
         */
        @Override
        public void remove() {
            throw new UnsupportedOperationException("Unable to remove from a BAMScheduleIterator");
        }

    }
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
    public final Map<SAMReaderID,GATKBAMFileSpan> fileSpans = new HashMap<SAMReaderID,GATKBAMFileSpan>();

    BAMScheduleEntry(final int start, final int stop) {
        this.start = start;
        this.stop = stop;
    }

    /**
     * Add a new file span to this schedule.
     * @param reader Reader associated with the span.
     * @param fileSpan Blocks to read in the given reader.
     */
    public void addFileSpan(final SAMReaderID reader, final GATKBAMFileSpan fileSpan) {
        fileSpans.put(reader,fileSpan);
    }

    /**
     * A naive merge operation.  Merge the fileSpans in other into this, blowing up if conflicts are
     * detected. Completely ignores merging start and stop.
     * @param other Other schedule entry to merging into this one.
     */
    public void mergeInto(final BAMScheduleEntry other) {
        final int thisSize = fileSpans.size();
        final int otherSize = other.fileSpans.size();
        fileSpans.putAll(other.fileSpans);
        if(fileSpans.size() != thisSize+otherSize)
            throw new ReviewedGATKException("Unable to handle overlaps when merging BAM schedule entries.");
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
