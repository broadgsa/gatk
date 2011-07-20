package org.broadinstitute.sting.alignment.reference.packing;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

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
    private final FileInputStream targetInputStream;

    /**
     * Channel source for packed bases.
     */
    private final FileChannel targetInputChannel;

    /**
     * A fixed-size buffer for word-packed data.
     */
    private final ByteOrder byteOrder;

    /**
     * How many bases are in a given packed word.
     */
    private final int basesPerPackedWord = PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BASE;

    /**
     * How many bytes in an integer?
     */
    private final int bytesPerInteger = PackUtils.bitsInType(Integer.class)/PackUtils.BITS_PER_BYTE;


    public BasePackedInputStream( Class<T> type, File inputFile, ByteOrder byteOrder ) throws FileNotFoundException {
        this(type,new FileInputStream(inputFile),byteOrder);
    }

    public BasePackedInputStream( Class<T> type, FileInputStream inputStream, ByteOrder byteOrder ) {
        if( type != Integer.class )
            throw new ReviewedStingException("Only bases packed into 32-bit words are currently supported by this input stream.  Type specified: " + type.getName());
        this.type = type;
        this.targetInputStream = inputStream;
        this.targetInputChannel = inputStream.getChannel();
        this.byteOrder = byteOrder;
    }

    /**
     * Read the entire contents of the input stream.
     * @param bwt array into which bases should be read.
     * @throws IOException if an I/O error occurs.
     */
    public void read(byte[] bwt) throws IOException {
        read(bwt,0,bwt.length);
    }

    /**
     * Read the next <code>length</code> bases into the bwt array, starting at the given offset.
     * @param bwt array holding the given data.
     * @param offset target position in the bases array into which bytes should be written.
     * @param length number of bases to read from the stream.
     * @throws IOException if an I/O error occurs.
     */
    public void read(byte[] bwt, int offset, int length) throws IOException {
        int bufferWidth = ((bwt.length+basesPerPackedWord-1)/basesPerPackedWord)*bytesPerInteger;
        ByteBuffer buffer = ByteBuffer.allocate(bufferWidth).order(byteOrder);
        targetInputChannel.read(buffer);
        targetInputChannel.position(targetInputChannel.position()+buffer.remaining());
        buffer.flip();

        int packedWord = 0;
        int i = 0;
        while(i < length) {
            if(i % basesPerPackedWord == 0) packedWord = buffer.getInt();
            int position = basesPerPackedWord - i%basesPerPackedWord - 1;
            bwt[offset+i++] = PackUtils.unpackBase((byte)((packedWord >> position*PackUtils.BITS_PER_BASE) & 0x3));            
        }
    }
}
