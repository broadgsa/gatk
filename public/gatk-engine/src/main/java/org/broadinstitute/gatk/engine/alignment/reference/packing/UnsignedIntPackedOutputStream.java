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

package org.broadinstitute.gatk.engine.alignment.reference.packing;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Writes an list of integers to the output file.
 *
 * @author mhanna
 * @version 0.1
 */
public class UnsignedIntPackedOutputStream {
    /**
     * Ultimate target for the occurrence array.
     */
    private final OutputStream targetOutputStream;

    /**
     * A fixed-size buffer for int-packed data.
     */
    private final ByteBuffer buffer;

    /**
     * Create a new PackedIntOutputStream, writing to the given target file.
     * @param outputFile target output file.
     * @param byteOrder Endianness to use when writing a list of integers.
     * @throws IOException if an I/O error occurs.
     */
    public UnsignedIntPackedOutputStream(File outputFile, ByteOrder byteOrder) throws IOException {
        this(new FileOutputStream(outputFile),byteOrder);
    }

    /**
     * Write packed ints to the given OutputStream.
     * @param outputStream Output stream to which to write packed ints.
     * @param byteOrder Endianness to use when writing a list of integers.
     */
    public UnsignedIntPackedOutputStream(OutputStream outputStream, ByteOrder byteOrder) {
        this.targetOutputStream = outputStream;
        buffer = ByteBuffer.allocate(PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BYTE).order(byteOrder);
    }

    /**
     * Write the data to the output stream.
     * @param datum datum to write. 
     * @throws IOException if an I/O error occurs.
     */
    public void write( long datum ) throws IOException {
        buffer.rewind();
        buffer.putInt((int)datum);
        targetOutputStream.write(buffer.array());
    }

    /**
     * Write the data to the output stream.
     * @param data data to write.  occurrences.length must match alphabet size.
     * @throws IOException if an I/O error occurs.
     */
    public void write( long[] data ) throws IOException {
        for(long datum: data)
            write(datum);
    }

    /**
     * Write the given chunk of data to the input stream.
     * @param data data to write.
     * @param offset position at which to start.
     * @param length number of ints to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( long[] data, int offset, int length ) throws IOException {
        for( int i = offset; i < offset+length; i++ )
            write(data[i]);
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
        targetOutputStream.close();
    }

}
