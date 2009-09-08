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

import org.broadinstitute.sting.utils.StingException;

import java.io.*;

/**
 * Write packed bases to an output stream.  Pack each base into 2 bits.
 *
 * @author mhanna
 * @version 0.1
 */
public class BytePackedOutputStream {
    /**
     * How many possible bases can be encoded?
     */
    public static final int ALPHABET_SIZE = 4;

    /**
     * Ultimate target for the packed bases.
     */
    private final OutputStream targetOutputStream;

    /**
     * The next byte to write to the output stream.  Will be added
     * to the output stream when enough bases are accumulated, or when
     * the file is closed.
     */
    private byte packedBases;

    /**
     * Where will the next base be embedded into packedBases?
     */
    private int positionInPack = 0;

    public BytePackedOutputStream( File outputFile ) throws FileNotFoundException {
        this(new BufferedOutputStream(new FileOutputStream(outputFile)));
    }

    /**
     * Write packed bases to the given output stream.
     * @param outputStream Output stream to which to write packed bases.
     */
    public BytePackedOutputStream( OutputStream outputStream ) {
        this.targetOutputStream = outputStream;
    }

    /**
     * Write a given base to the output stream.
     * @param base Base to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( byte base ) throws IOException {
        packedBases |= (getPackedRepresentation(base) << 2*(ALPHABET_SIZE-positionInPack-1));

        // Increment the packed counter.  If all possible bases have been squeezed into this byte, write it out.
        positionInPack = ++positionInPack % ALPHABET_SIZE;
        if( positionInPack == 0 ) {
            targetOutputStream.write(packedBases);
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
            targetOutputStream.write(packedBases);
            targetOutputStream.write(positionInPack);
        }
        else
            targetOutputStream.write(ALPHABET_SIZE);

        targetOutputStream.close();
    }

    /**
     * Gets the two-bit representation of a base.  A=00b, C=01b, G=10b, T=11b.
     * @param base ASCII value for the base to pack.
     * @return A byte from 0-3 indicating the base's packed value.
     */
    public static byte getPackedRepresentation(byte base) {
        switch( base ) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                throw new StingException("Unknown base type: " + base);
        }
    }

    public static byte decodePackedRepresentation(byte pack) {
        switch( pack ) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                throw new StingException("Unknown pack type: " + pack);
        }
    }

}
