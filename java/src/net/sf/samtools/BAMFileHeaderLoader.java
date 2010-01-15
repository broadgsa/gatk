package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.StringLineReader;

import java.io.File;
import java.io.IOException;
import java.io.DataInputStream;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

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

    public static final Chunk preambleLocation = new Chunk(0<<16 | 0, 0<<16 | 3);

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

        // directly copied from BAMFileReader...
        final int sequenceCount = binaryCodec.readInt();
        if (header.getSequenceDictionary().size() > 0) {
            // It is allowed to have binary sequences but no text sequences, so only validate if both are present
            if (sequenceCount != header.getSequenceDictionary().size()) {
                throw new SAMFormatException("Number of sequences in text header (" +
                        header.getSequenceDictionary().size() +
                        ") != number of sequences in binary header (" + sequenceCount + ") for file " + file);
            }
            for (int i = 0; i < sequenceCount; i++) {
                final SAMSequenceRecord binarySequenceRecord = readSequenceRecord(binaryCodec,file);
                final SAMSequenceRecord sequenceRecord = header.getSequence(i);
                if (!sequenceRecord.getSequenceName().equals(binarySequenceRecord.getSequenceName())) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different names in file " +
                            binaryCodec);
                }
                if (sequenceRecord.getSequenceLength() != binarySequenceRecord.getSequenceLength()) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different lengths in file " +
                            binaryCodec);
                }
            }
        } else {
            // If only binary sequences are present, copy them into mFileHeader
            final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>(sequenceCount);
            for (int i = 0; i < sequenceCount; i++) {
                sequences.add(readSequenceRecord(binaryCodec,file));
            }
            header.setSequenceDictionary(new SAMSequenceDictionary(sequences));
        }
        inputStream.close();

        return new BAMFileHeaderLoader(header,new Chunk(buffer.length,inputStream.getFilePointer()-1));
    }

    /**
     * Reads a single binary sequence record from the file or stream
     * @param binaryCodec stream to read from.
     * @param file  Note that this is used only for reporting errors.
     * @return an individual sequence record.
     */
    private static SAMSequenceRecord readSequenceRecord(final BinaryCodec binaryCodec, final File file) {
        final int nameLength = binaryCodec.readInt();
        if (nameLength <= 1) {
            throw new SAMFormatException("Invalid BAM file header: missing sequence name in file " + file.getAbsolutePath());
        }
        final String sequenceName = binaryCodec.readString(nameLength - 1);
        // Skip the null terminator
        binaryCodec.readByte();
        final int sequenceLength = binaryCodec.readInt();
        return new SAMSequenceRecord(sequenceName, sequenceLength);
    }
}
