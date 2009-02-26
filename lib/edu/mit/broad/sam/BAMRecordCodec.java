/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.RuntimeEOFException;
import edu.mit.broad.sam.util.SortingCollection;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.Map;

public class BAMRecordCodec implements SortingCollection.Codec<SAMRecord> {
    private final BinaryCigarCodec cigarCodec = new BinaryCigarCodec();
    private final SAMFileHeader header;
    private OutputStream os;
    private InputStream is;
    private BinaryCodec binaryCodec;
    private BinaryTagCodec binaryTagCodec;

    public BAMRecordCodec(final SAMFileHeader header) {
        this.header = header;
    }

    public BAMRecordCodec clone() {
        BAMRecordCodec other = new BAMRecordCodec(this.header);
        return other;
    }


    /** Sets the output stream that records will be written to. */
    public void setOutputStream(final OutputStream os) {
        this.os = os;
        this.binaryCodec    = new BinaryCodec(this.os);
        this.binaryTagCodec = new BinaryTagCodec(this.binaryCodec);
    }

    /** Sets the input stream that records will be read from. */
    public void setInputStream(final InputStream is) {
        this.is = is;
        this.binaryCodec    = new BinaryCodec(this.is);
        this.binaryTagCodec = new BinaryTagCodec(this.binaryCodec);
    }

    /**
     * Write object to OutputStream.
     * The SAMRecord must have a header set into it so reference indices can be resolved.
     *
     * @param alignment what to write
     */
    public void encode(final SAMRecord alignment) {
        // Compute block size, as it is the first element of the file representation of SAMRecord
        final int readLength = alignment.getReadLength();

        final int cigarLength = alignment.getCigarLength();

        int blockSize = BAMFileConstants.FIXED_BLOCK_SIZE + alignment.getReadNameLength() + 1  + // null terminated
                        cigarLength * 4 +
                        (readLength + 1) / 2 + // 2 bases per byte
                        readLength;

        final int attributesSize = alignment.getAttributesBinarySize();
        if (attributesSize != -1) {
            blockSize += attributesSize;
        } else {
            if (alignment.getAttributes() != null) {
                for (final Map.Entry<String, Object> attribute : alignment.getAttributes()) {
                    blockSize += (BinaryTagCodec.getTagSize(attribute.getValue()));
                }
            }
        }

        int indexBin = 0;
        if (alignment.getReferenceIndex(header) >= 0) {
            if (alignment.getIndexingBin() != null) {
                indexBin = alignment.getIndexingBin();
            } else {
                indexBin = SAMUtils.reg2bin(alignment.getAlignmentStart() - 1,
                        alignment.getAlignmentEnd() - 1);
            }
        }

        // Blurt out the elements
        this.binaryCodec.writeInt(blockSize);
        this.binaryCodec.writeInt(alignment.getReferenceIndex(header));
        // 0-based!!
        this.binaryCodec.writeInt(alignment.getAlignmentStart() - 1);
        this.binaryCodec.writeUByte((short)(alignment.getReadNameLength() + 1));
        this.binaryCodec.writeUByte((short)alignment.getMappingQuality());
        this.binaryCodec.writeUShort(indexBin);
        this.binaryCodec.writeUShort(cigarLength);
        this.binaryCodec.writeUShort(alignment.getFlags());
        this.binaryCodec.writeInt(alignment.getReadLength());
        this.binaryCodec.writeInt(alignment.getMateReferenceIndex(header));
        this.binaryCodec.writeInt(alignment.getMateAlignmentStart() - 1);
        this.binaryCodec.writeInt(alignment.getInferredInsertSize());
        final byte[] variableLengthBinaryBlock = alignment.getVariableBinaryRepresentation();
        if (variableLengthBinaryBlock != null) {
            this.binaryCodec.writeBytes(variableLengthBinaryBlock);
        } else {
            this.binaryCodec.writeString(alignment.getReadName(), false, true);
            final int[] binaryCigar = cigarCodec.encode(alignment.getCigar());
            for (final int cigarElement : binaryCigar) {
                // Assumption that this will fit into an integer, despite the fact
                // that it is specced as a uint.
                this.binaryCodec.writeInt(cigarElement);
            }
            this.binaryCodec.writeBytes(SAMUtils.bytesToCompressedBases(alignment.getReadBases()));
            this.binaryCodec.writeBytes(alignment.getBaseQualities());
            if (alignment.getAttributes() != null) {
                for (final Map.Entry<String, Object> attribute : alignment.getAttributes()) {
                    this.binaryTagCodec.writeTag(attribute.getKey(), attribute.getValue());
                }
            }
        }
    }

    /**
     * Read the next record from the input stream and convert into a java object.
     *
     * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
     *         a record.
     */
    public SAMRecord decode() {
        int recordLength = 0;
        try {
            recordLength = this.binaryCodec.readInt();
        }
        catch (RuntimeEOFException e) {
            return null;
        }

        if (recordLength < BAMFileConstants.FIXED_BLOCK_SIZE ||
                recordLength > BAMFileConstants.MAXIMUM_RECORD_LENGTH) {
            throw new SAMFormatException("Invalid record length: " + recordLength);
        }
        
        final int referenceID = this.binaryCodec.readInt();
        final int coordinate = this.binaryCodec.readInt() + 1;
        final short readNameLength = this.binaryCodec.readUByte();
        final short mappingQuality = this.binaryCodec.readUByte();
        final int bin = this.binaryCodec.readUShort();
        final int cigarLen = this.binaryCodec.readUShort();
        final int flags = this.binaryCodec.readUShort();
        final int readLen = this.binaryCodec.readInt();
        final int mateReferenceID = this.binaryCodec.readInt();
        final int mateCoordinate = this.binaryCodec.readInt() + 1;
        final int insertSize = this.binaryCodec.readInt();
        final byte[] restOfRecord = new byte[recordLength - BAMFileConstants.FIXED_BLOCK_SIZE];
        this.binaryCodec.readBytes(restOfRecord);
        final BAMRecord ret = new BAMRecord(header, referenceID, coordinate, readNameLength, mappingQuality,
                bin, cigarLen, flags, readLen, mateReferenceID, mateCoordinate, insertSize, restOfRecord);
        ret.setHeader(header); 
        return ret;
    }
}
