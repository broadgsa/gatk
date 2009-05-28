package org.broadinstitute.sting.utils.glf;

import net.sf.samtools.util.BinaryCodec;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: May 13, 2009
 * Time: 3:36:18 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 13, 2009
 * <p/>
 * Class GLFRecord
 * <p/>
 * The base record that's stored in the GLF format.
 */
public class GLFRecord {
    // our record list
    private final List<RecordType> records = new ArrayList<RecordType>();

    public static final byte[] glfMagic = {'G','L','F','\3'};
    private String headerText = "";
    private String referenceSequenceName = "";
    private long referenceSequenceLength = 0;

    private int currentOffset = -1;
    /**
     * The public constructor for creating a GLF object
     * @param headerText the header text (currently unclear what the contents are)
     * @param referenceSequenceName the reference sequence name
     */
    public GLFRecord(String headerText, String referenceSequenceName, int referenceSequenceLength) {
        this.headerText = headerText;
        this.referenceSequenceName = referenceSequenceName;
        this.referenceSequenceLength = referenceSequenceLength;
    }

    public void addSNPCall(int genomicLoc, long read_depth, int rmsMapQ, LikelihoodObject lhValues) {
        if (currentOffset >= genomicLoc) {
            throw new IllegalArgumentException("The location supplied is less then a previous location");
        }

        // make sure the read depth isn't too large
        if (read_depth < 0 || read_depth > 0x00FFFFFF) {
            throw new IllegalArgumentException("The read depth is too large; must lie in the range 0 to 0x00ffffff");
        }

        // check that the rmsSquare is greater then 0, and will fit in a uint8
        if (rmsMapQ > 0x000000FF || rmsMapQ < 0) {
            throw new IllegalArgumentException("rms of mapping quality is too large; must lie in the range 0 to 0x000000ff");
        }

        if (lhValues == null) {
            throw new IllegalArgumentException("likelihood object cannot be null");
        }

        SinglePointCall call = new SinglePointCall(genomicLoc - currentOffset,
                                                   read_depth,
                                                   rmsMapQ,
                                                   lhValues.toByteArray(),
                                                   (short)lhValues.getMinimumValue());

    }

    /**
     * Write out the record to a binary codec object
     *
     * @param out
     */
    public void write(BinaryCodec out) {
        out.writeBytes(glfMagic);
        out.writeString(headerText,true,true);
        out.writeString(referenceSequenceName,true,true);
        out.writeUInt(referenceSequenceLength);
        for (RecordType rec: records) {
            out.writeUByte(rec.getRecordType().getFieldValue());
            rec.write(out);
        }
        out.writeByte((byte)0);
    }

}


interface RecordType {
    
    enum RECORD_TYPE {
        VARIABLE((short)2),
        SINGLE((short)1);
        private final short fieldValue;   // in kilograms
        RECORD_TYPE(short value) {
            fieldValue = value;
        }
        public short getFieldValue() {
            return fieldValue;
        }
    };

    /**
     * write the record out to a binary codec
     * @param out
     */
    public void write(BinaryCodec out);

    /**
     * get the record type
     * @return the record type
     */
    public RECORD_TYPE getRecordType();

    /**
     * 
     * @return
     */
    public int getByteSize();
}

// the second record type
class VariableLengthCall implements RecordType {
    public int offset = 0;
    public int min_depth = 0;
    public byte rmsMapQ = 0;
    public byte lkHom1 = 0;
    public byte lkHom2 = 0;
    public byte lkHet = 0;
    public short indelLen1 = 0;
    public short indelLen2 = 0;
    public final byte indelSeq1[];
    public final byte indelSeq2[];

    // our size, we're immutable (at least the size is)
    private final int size; // in bytes

    /**
     * the BinaryCodec constructor
     *
     * @param in the binary codec to get data from
     */
    VariableLengthCall(BinaryCodec in) {
        offset = in.readInt();
        min_depth = in.readInt();
        rmsMapQ = in.readByte();
        lkHom1 = in.readByte();
        lkHom2 = in.readByte();
        lkHet = in.readByte();
        indelLen1 = in.readShort();
        indelLen2 = in.readShort();
        indelSeq1 = new byte[indelLen1];
        indelSeq2 = new byte[indelLen2];
        this.size = 16 + indelLen1 + indelLen2; 
    }

    /**
     * Write out the record to a binary codec
     *
     * @param out
     */
    public void write(BinaryCodec out) {
        out.writeInt(offset);
        out.writeInt(min_depth);
        out.writeByte(rmsMapQ);
        out.writeByte(lkHom1);
        out.writeByte(lkHom2);
        out.writeByte(lkHet);
        out.writeShort(indelLen1);
        out.writeShort(indelLen2);
        out.writeBytes(indelSeq1);
        out.writeBytes(indelSeq2);
    }

    public RECORD_TYPE getRecordType() {
        return RECORD_TYPE.VARIABLE;
    }

    /** @return  */
    public int getByteSize() {
        return size;
    }
}


// the first record type
class SinglePointCall implements RecordType {
    // our likelyhood array size
    public static final int LIKELYHOOD_SIZE = 10;

    // class fields
    public int offset = 0;
    public Long min_depth = 0L;
    public int rmsMapQ = 0;
    public final short lk[] = new short[LIKELYHOOD_SIZE];
    public short minimumLikelihood = 0;

    // our size, we're immutable (the size at least)
    private final int byteSize; // in bytes

    /**
     * The parameter constructor
     *
     * @param offset
     * @param min_depth
     * @param rmsMapQ
     * @param lk
     */
    SinglePointCall(int offset, long min_depth, int rmsMapQ, short[] lk, short minimumLikelihood) {
        if (lk.length != LIKELYHOOD_SIZE) {
            throw new IllegalArgumentException("SinglePointCall: passed in likelyhood array size != LIKELYHOOD_SIZE");
        }
        this.offset = offset;
        this.min_depth = min_depth;
        this.rmsMapQ = rmsMapQ;
        this.minimumLikelihood = minimumLikelihood;
        System.arraycopy(lk, 0, this.lk, 0, LIKELYHOOD_SIZE);
        byteSize = 9 + lk.length;
    }

    /**
     * the BinaryCodec constructor
     *
     * @param in the binary codec to get data from
     */
    SinglePointCall(BinaryCodec in) {
        offset = in.readInt();
        min_depth = Long.valueOf(in.readInt());
        rmsMapQ = in.readByte();
        for (int x = 0; x < LIKELYHOOD_SIZE; x++) {
            lk[x] = in.readUByte();
        }
        byteSize = 9 + lk.length;
    }

    /**
     * Write out the record to a binary codec
     *
     * @param out
     */
    public void write(BinaryCodec out) {
        out.writeInt(offset);
        out.writeInt(min_depth.intValue());
        out.writeByte(rmsMapQ);
        for (int x = 0; x < LIKELYHOOD_SIZE; x++) {
             out.writeUByte(lk[x]);
        }
    }

    public RECORD_TYPE getRecordType() {
        return RECORD_TYPE.SINGLE;
    }

    /** @return  */
    public int getByteSize() {
        return byteSize;
    }

}