package org.broadinstitute.sting.bwa;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.ByteBuffer;
import java.util.List;
import java.util.ArrayList;

/**
 * Reads a word-packed version of the input stream.
 *
 * @author mhanna
 * @version 0.1
 */
public class WordPackedInputStream {

    /**
     * Ultimate source for packed bases.
     */
    private final InputStream targetInputStream;

    /**
     * A fixed-size buffer for word-packed data.
     */
    private final ByteBuffer buffer;

    public WordPackedInputStream( File inputFile, ByteOrder byteOrder ) throws FileNotFoundException {
        this.targetInputStream = new BufferedInputStream(new FileInputStream(inputFile));
        this.buffer = ByteBuffer.allocate(WordPackedOutputStream.BASES_PER_WORD/BytePackedOutputStream.ALPHABET_SIZE).order(byteOrder);
    }

    /**
     * Read the entire contents of the input stream.
     * @throws IOException if an I/O error occurs.
     */
    public byte[] read() throws IOException {
        // Skip over header info.
        for( int i = 0; i < 5; i++ ) {
            targetInputStream.read(buffer.array());
            System.out.println("Skipping over: " + buffer.getInt());
            buffer.rewind();
        }

        List<Byte> bwtList = new ArrayList<Byte>();
        while(targetInputStream.read(buffer.array()) > 0) {
            int packedWord = buffer.getInt();
            for( int i = WordPackedOutputStream.BASES_PER_WORD-1; i >= 0; i-- ) {
                byte packedByte = (byte)((packedWord >> i*2) & 0x3);
                bwtList.add(BytePackedOutputStream.decodePackedRepresentation(packedByte));
            }
            buffer.rewind();
        }

        byte[] bwt = new byte[bwtList.size()];
        for(int i = 0; i < bwtList.size(); i++)
            bwt[i] = bwtList.get(i);

        return bwt;
    }

}
