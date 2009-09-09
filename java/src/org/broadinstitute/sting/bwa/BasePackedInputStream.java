package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.ByteBuffer;

/**
 * Reads a packed version of the input stream.
 *
 * @author mhanna
 * @version 0.1
 */
public class BasePackedInputStream<T> {
    /**
     * Type of object to unpack.
     */
    private final Class<T> type;

    /**
     * Ultimate source for packed bases.
     */
    private final InputStream targetInputStream;

    /**
     * A fixed-size buffer for word-packed data.
     */
    private final ByteBuffer buffer;

    public BasePackedInputStream( Class<T> type, File inputFile, ByteOrder byteOrder ) throws FileNotFoundException {
        this(type,new BufferedInputStream(new FileInputStream(inputFile)),byteOrder);
    }

    public BasePackedInputStream( Class<T> type, InputStream inputStream, ByteOrder byteOrder ) {
        if( type != Integer.class )
            throw new StingException("Only bases packed into 32-bit words are currently supported by this input stream.  Type specified: " + type.getName());

        this.targetInputStream = inputStream;
        this.type = type;
        this.buffer = ByteBuffer.allocate(PackUtils.bitsInType(type)/PackUtils.BITS_PER_BYTE).order(byteOrder);
    }

    /**
     * Read the entire contents of the input stream.
     * @param length number of bases to read from the stream.
     * @return a byte array of the given length.
     * @throws IOException if an I/O error occurs.
     */
    public byte[] read( int length ) throws IOException {
        byte[] bwt = new byte[length];
        int packedWord = 0;

        final int basesPerEntry = PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BASE;
        for( int i = 0; i < length; i++ ) {
            if( i % basesPerEntry == 0 ) {
                buffer.rewind();
                targetInputStream.read(buffer.array());
                packedWord = buffer.getInt();
            }

            int position = basesPerEntry - i % basesPerEntry - 1;
            bwt[i] = PackUtils.unpackBase((byte)((packedWord >> position*PackUtils.BITS_PER_BASE) & 0x3));
        }

        return bwt;
    }

}
