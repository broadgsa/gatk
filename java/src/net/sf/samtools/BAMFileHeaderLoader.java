package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.StringLineReader;

import java.io.File;
import java.io.IOException;
import java.io.DataInputStream;
import java.util.Arrays;

/**
 * Loads a BAM file header from an file, optionally providing its position
 * within the file.
 *
 * @author mhanna
 * @version 0.1
 */
public class BAMFileHeaderLoader {
    /**
     * The contents of the BAM file header.
     */
    private final SAMFileHeader header;

    /**
     * Location of the header within the BAM file.
     */
    private final Chunk location;

    /**
     * Load the header from the given file.
     * @param header the parsed haeder for the BAM file.
     * @param location the location of the header (start and stop) within the BAM.
     */
    private BAMFileHeaderLoader(SAMFileHeader header, Chunk location) {
        this.header = header;
        this.location = location;
    }

    /**
     * Gets the header for the given BAM file.
     * @return The header for this BAM file.
     */
    public SAMFileHeader getHeader() {
        return header;
    }

    /**
     * Gets the location of the header within the given BAM file, in chunk format.
     * @return the location of the header, in chunk coordinates.
     */
    public Chunk getLocation() {
        return location;
    }

    public static BAMFileHeaderLoader load(File file) throws IOException {
        BlockCompressedInputStream inputStream = new BlockCompressedInputStream(file);
        BinaryCodec binaryCodec = new BinaryCodec(new DataInputStream(inputStream));

        final byte[] buffer = new byte[4];
        binaryCodec.readBytes(buffer);
        if (!Arrays.equals(buffer, BAMFileConstants.BAM_MAGIC)) {
            throw new IOException("Invalid BAM file header");
        }

        final int headerTextLength = binaryCodec.readInt();
        final String textHeader = binaryCodec.readString(headerTextLength);
        final SAMTextHeaderCodec headerCodec = new SAMTextHeaderCodec();
        headerCodec.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        SAMFileHeader header = headerCodec.decode(new StringLineReader(textHeader),file.getAbsolutePath());

        inputStream.close();

        return new BAMFileHeaderLoader(header,new Chunk(buffer.length,inputStream.getFilePointer()));
    }

}
