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

import net.sf.samtools.GATKBin;
import net.sf.samtools.GATKChunk;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Creates an 'index of the index' for a particular reference sequence
 * within the BAM file, for easier whole-BAM-file traversal.
 * file.
 */
public class BAMIndexBinIterator {
    /**
     * The source of index data.
     */
    private final GATKBAMIndex index;

    /**
     * The file storing the index data.
     */
    private final File indexFile;

    /**
     * File storing index metadata.
     */
    private final File metaIndexFile;

    /**
     * Reference sequence that temporary file is based on.
     */
    private final int referenceSequence;

    /**
     * Size of a long in bytes.
     */
    private static final int LONG_SIZE_IN_BYTES = Long.SIZE / 8;

    public BAMIndexBinIterator(final GATKBAMIndex index, final File indexFile, final int referenceSequence) {
        this.index = index;
        this.indexFile = indexFile;
        this.referenceSequence = referenceSequence;

        index.seek(4);

        final int sequenceCount = index.readInteger();

        if (referenceSequence >= sequenceCount)
            throw new ReviewedStingException(String.format("Reference sequence past end of genome; reference sequence = %d, sequence count = %d",referenceSequence,sequenceCount));

        index.skipToSequence(referenceSequence);

        int binCount = index.readInteger();

        try {
            metaIndexFile = File.createTempFile("bammetaindex."+referenceSequence,null);
            metaIndexFile.deleteOnExit();

            FileOutputStream metaIndex = new FileOutputStream(metaIndexFile);
            FileChannel metaIndexChannel = metaIndex.getChannel();

            // zero out the contents of the file.  Arrays of primitives in java are always zeroed out by default.
            byte[] emptyContents = new byte[GATKBAMIndex.MAX_BINS*(Long.SIZE/8)]; // byte array is zeroed out by default.
            metaIndexChannel.write(ByteBuffer.wrap(emptyContents));

            ByteBuffer binPositionBuffer = ByteBuffer.allocate(emptyContents.length);
            binPositionBuffer.order(ByteOrder.LITTLE_ENDIAN);            

            for (int binNumber = 0; binNumber < binCount; binNumber++) {
                long position = index.position();

                final int indexBin = index.readInteger();

                metaIndexChannel.position(indexBin*LONG_SIZE_IN_BYTES);
                binPositionBuffer.putLong(position);
                binPositionBuffer.flip();

                System.out.printf("Writing bin number %d to position %d: coordinate = %d%n",indexBin,indexBin*Long.SIZE*8,position);

                metaIndexChannel.write(binPositionBuffer);
                binPositionBuffer.flip();

                final int nChunks = index.readInteger();
                index.skipBytes(16 * nChunks);
            }

            metaIndexChannel.close();
            metaIndex.close();
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to write BAM metaindex",ex);
        }
    }

    public void close() {
        metaIndexFile.delete();    
    }

    public CloseableIterator<GATKBin> getIteratorOverLevel(int level) {
        return new LevelIterator(level);
    }

    private class LevelIterator implements CloseableIterator<GATKBin> {
        /**
         * The raw BAM index file with unordered bins.
         */
        private final FileInputStream indexInputStream;

        /**
         * The index metafile, with pointers to ordered bins.
         */
        private final FileInputStream metaIndexInputStream;

        /**
         * The first and last bins in the level.
         */
        private final int firstBinNumber, lastBinNumber;

        /**
         * The current bin in the index.
         */
        private int currentBinNumber;

        /**
         * Position of the most recent chunk data.
         */
        private GATKBin nextBin = null;

        public LevelIterator(final int level) {
            try {
                indexInputStream = new FileInputStream(indexFile);
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to open index file for reading");
            }
            try {
                metaIndexInputStream = new FileInputStream(metaIndexFile);
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to open index metafile for reading");
            }

            firstBinNumber = GATKBAMIndex.getFirstBinInLevel(level);
            lastBinNumber = firstBinNumber + index.getLevelSize(level) - 1;

            currentBinNumber = firstBinNumber - 1;
            advance();
        }

        public void close() {
            try {
                indexInputStream.close();
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to close index file.");
            }
            try {
                metaIndexInputStream.close();
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to close index metafile");
            }
        }

        public boolean hasNext() {
            return nextBin != null;
        }

        public GATKBin next() {
            if(!hasNext())
                throw new NoSuchElementException("Out of elements in BAMIndexBinIterator");
            GATKBin currentPosition = nextBin;
            advance();
            return currentPosition;
        }

        public void remove() { throw new UnsupportedOperationException("Cannot remove from a LevelIterator"); }

        private void advance() {
            ByteBuffer indexFilePositionDecoder = ByteBuffer.allocate(LONG_SIZE_IN_BYTES*2);
            indexFilePositionDecoder.order(ByteOrder.LITTLE_ENDIAN);

            nextBin = null;
            try {
                indexFilePositionDecoder.limit(LONG_SIZE_IN_BYTES);
                while(nextBin == null && currentBinNumber < lastBinNumber) {
                    currentBinNumber++;
                    metaIndexInputStream.getChannel().position(currentBinNumber*LONG_SIZE_IN_BYTES);
                    metaIndexInputStream.getChannel().read(indexFilePositionDecoder);
                    indexFilePositionDecoder.flip();
                    long currentPosition = indexFilePositionDecoder.getLong();
                    indexFilePositionDecoder.flip();

                    if(currentPosition != 0) {
                        indexInputStream.getChannel().position(currentPosition);
                        indexInputStream.getChannel().read(indexFilePositionDecoder);

                        indexFilePositionDecoder.flip();
                        int binNumber = indexFilePositionDecoder.getInt();
                        if(binNumber != currentBinNumber)
                            throw new ReviewedStingException("Index file and metaindex file are out of sync.");

                        int nChunks = indexFilePositionDecoder.getInt();
                        GATKChunk[] chunks = new GATKChunk[nChunks];

                        indexFilePositionDecoder.limit(LONG_SIZE_IN_BYTES*2);
                        indexFilePositionDecoder.clear();

                        for (int ci = 0; ci < nChunks; ci++) {
                            indexInputStream.getChannel().read(indexFilePositionDecoder);

                            indexFilePositionDecoder.flip();
                            final long chunkBegin = indexFilePositionDecoder.getLong();
                            final long chunkEnd = indexFilePositionDecoder.getLong();
                            chunks[ci] = new GATKChunk(chunkBegin, chunkEnd);

                            indexFilePositionDecoder.flip();                            
                        }

                        nextBin = new GATKBin(referenceSequence,binNumber);
                        nextBin.setChunkList(chunks);
                    }
                }
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to close index metafile");
            }

        }
    }
}
