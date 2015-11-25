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
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Read a set of integers packed into 
 *
 * @author mhanna
 * @version 0.1
 */
public class UnsignedIntPackedInputStream {
    /**
     * Ultimate target for the occurrence array.
     */
    private final FileInputStream targetInputStream;

    /**
     * Target channel from which to pull file data.
     */
    private final FileChannel targetInputChannel;

    /**
     * The byte order in which integer input data appears.
     */
    private final ByteOrder byteOrder;

    /**
     * How many bytes are required to store an integer?
     */
    private final int bytesPerInteger = PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BYTE;

    /**
     * Create a new PackedIntInputStream, writing to the given target file.
     * @param inputFile target input file.
     * @param byteOrder Endianness to use when writing a list of integers.
     * @throws java.io.IOException if an I/O error occurs.
     */
    public UnsignedIntPackedInputStream(File inputFile, ByteOrder byteOrder) throws IOException {
        this(new FileInputStream(inputFile),byteOrder);
    }

    /**
     * Read  ints from the given InputStream.
     * @param inputStream Input stream from which to read ints.
     * @param byteOrder Endianness to use when writing a list of integers.
     */
    public UnsignedIntPackedInputStream(FileInputStream inputStream, ByteOrder byteOrder) {
        this.targetInputStream = inputStream;
        this.targetInputChannel = inputStream.getChannel();
        this.byteOrder = byteOrder;
    }

    /**
     * Read a datum from the input stream.
     * @return The next input datum in the stream.
     * @throws IOException if an I/O error occurs.
     */
    public long read() throws IOException {
        long[] data = new long[1];
        read(data);
        return data[0];
    }

    /**
     * Read the data from the input stream.
     * @param data placeholder for input data.
     * @throws IOException if an I/O error occurs.
     */
    public void read( long[] data ) throws IOException {
        read( data, 0, data.length );
    }

    /**
     * Read the data from the input stream, starting at the given offset.
     * @param data placeholder for input data.
     * @param offset place in the array to start reading in data.
     * @param length number of ints to read in. 
     * @throws IOException if an I/O error occurs.
     */
    public void read( long[] data, int offset, int length ) throws IOException {
        ByteBuffer readBuffer = ByteBuffer.allocate(bytesPerInteger*length).order(byteOrder);

        targetInputChannel.read(readBuffer,targetInputChannel.position());
        readBuffer.flip();
        targetInputChannel.position(targetInputChannel.position()+readBuffer.remaining());

        int i = 0;
        while(i < length)
            data[offset+i++] = readBuffer.getInt() & 0xFFFFFFFFL;
    }

    /**
     * Closes the given output stream.
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        targetInputStream.close();
    }
}
