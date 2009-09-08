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

package org.broadinstitute.sting.bwa;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Word-pack bases into the output stream.  Bytes are stored as
 * little-endian unsigned ints.
 *
 * @author mhanna
 * @version 0.1
 */
public class WordPackedOutputStream {
    /**
     * How many bases can be stored in the given word?
     */
    public static final int BASES_PER_WORD = 16;

    /**
     * Ultimate target for the packed bases.
     */
    private final OutputStream targetOutputStream;

    /**
     * The next byte to write to the output stream.  Will be added
     * to the output stream when enough bases are accumulated, or when
     * the file is closed.
     */
    private int packedBases;

    /**
     * Where will the next base be embedded into packedBases?
     */
    private int positionInPack = 0;

    /**
     * A fixed-size buffer for word-packed data.
     */
    private final ByteBuffer buffer;

    public WordPackedOutputStream( File outputFile, ByteOrder byteOrder ) throws FileNotFoundException {
        this(new BufferedOutputStream(new FileOutputStream(outputFile)),byteOrder);
    }

    /**
     * Write packed bases to the given output stream.
     * @param outputStream Output stream to which to write packed bases.
     * @param byteOrder Switch between big endian / little endian when reading / writing files.
     */
    public WordPackedOutputStream(OutputStream outputStream, ByteOrder byteOrder) {
        this.targetOutputStream = outputStream;
        this.buffer = ByteBuffer.allocate(BASES_PER_WORD/BytePackedOutputStream.ALPHABET_SIZE).order(byteOrder);
    }

    /**
     * Write a given base to the output stream.
     * @param base Base to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( byte base ) throws IOException {
        packedBases |= (BytePackedOutputStream.getPackedRepresentation(base) << 2*(BASES_PER_WORD-positionInPack-1));

        // Increment the packed counter.  If all possible bases have been squeezed into this byte, write it out.
        positionInPack = ++positionInPack % BASES_PER_WORD;
        if( positionInPack == 0 ) {
            buffer.rewind();
            buffer.putInt(packedBases);
            targetOutputStream.write(buffer.array());
            packedBases = 0;
        }
    }

    /**
     * Writes an array of bases to the target output stream.
     * @param bases List of bases to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( byte[] bases ) throws IOException {
        for(byte base: bases) write(base);
    }

    /**
     * Flush the contents of the OutputStream to disk.
     * @throws IOException if an I/O error occurs.
     */
    public void flush() throws IOException {
        targetOutputStream.flush();
    }

    /**
     * Closes the given output stream.
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        // Write (incomplete) block in file, and number of bases in that last byte.
        if( positionInPack > 0 ) {
            buffer.rewind();
            buffer.putInt(packedBases);
            targetOutputStream.write(buffer.array());
        }
        targetOutputStream.close();
    }

}

