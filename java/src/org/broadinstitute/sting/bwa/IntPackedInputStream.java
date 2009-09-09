package org.broadinstitute.sting.bwa;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Read a set of integers packed into 
 *
 * @author mhanna
 * @version 0.1
 */
public class IntPackedInputStream {
    /**
     * Ultimate target for the occurrence array.
     */
    private final InputStream targetInputStream;

    /**
     * The byte order in which integer input data appears.
     */
    private final ByteBuffer buffer;

    /**
     * Create a new PackedIntInputStream, writing to the given target file.
     * @param inputFile target input file.
     * @param byteOrder Endianness to use when writing a list of integers.
     * @throws java.io.IOException if an I/O error occurs.
     */
    public IntPackedInputStream(File inputFile, ByteOrder byteOrder) throws IOException {
        this(new FileInputStream(inputFile),byteOrder);
    }

    /**
     * Read  ints from the given InputStream.
     * @param inputStream Input stream from which to read ints.
     * @param byteOrder Endianness to use when writing a list of integers.
     * @throws IOException if an I/O error occurs.
     */
    public IntPackedInputStream(InputStream inputStream, ByteOrder byteOrder) throws IOException {
        this.targetInputStream = inputStream;
        this.buffer = ByteBuffer.allocate(PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BYTE).order(byteOrder);        
    }

    /**
     * Read a datum from the input stream.
     * @return The next input datum in the stream.
     * @throws IOException if an I/O error occurs.
     */
    public int read() throws IOException {
        int[] data = new int[1];
        read(data);
        return data[0];
    }

    /**
     * Read the data from the input stream.
     * @param data placeholder for input data.
     * @throws IOException if an I/O error occurs.
     */
    public void read( int[] data ) throws IOException {
        for(int i = 0; i < data.length; i++) {
            targetInputStream.read(buffer.array());
            data[i] = buffer.getInt();
            buffer.rewind();
        }
    }

    /**
     * Closes the given output stream.
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        targetInputStream.close();
    }
}
