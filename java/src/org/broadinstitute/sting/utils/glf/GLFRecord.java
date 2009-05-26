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

    /**
     * The public constructor for creating a GLF object
     * @param headerText the header text (currently unclear what the contents are)
     * @param referenceSequenceName the reference sequence name
     */
    public GLFRecord(String headerText, String referenceSequenceName) {
        this.headerText = headerText;
        this.referenceSequenceName = referenceSequenceName;
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

    public void write(BinaryCodec out);

    public RECORD_TYPE getRecordType();

}

// the second record type
class VariableLengthGenotype implements RecordType {
    public int offset;
    public int min_depth;
    public byte rmsMapQ;
    public byte lkHom1;
    public byte lkHom2;
    public byte lkHet;
    public short indelLen1;
    public short indelLen2;
    public final byte indelSeq1[];
    public final byte indelSeq2[];

    // our size, we're immutable (at least the size is)
    private final int size; // in bytes

    /**
     * the BinaryCodec constructor
     *
     * @param in the binary codec to get data from
     */
    VariableLengthGenotype(BinaryCodec in) {
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
}


// the first record type
class SinglePointGenotype implements RecordType {
    // our likelyhood array size
    public static final int LIKELYHOOD_SIZE = 10;

    // class fields
    public int offset;
    public int min_depth;
    public byte rmsMapQ;
    public final byte lk[] = new byte[LIKELYHOOD_SIZE];

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
    SinglePointGenotype(int offset, int min_depth, byte rmsMapQ, byte[] lk) {
        if (lk.length != LIKELYHOOD_SIZE) {
            throw new IllegalArgumentException("SinglePointGenotype: passed in likelyhood array size != LIKELYHOOD_SIZE");
        }
        this.offset = offset;
        this.min_depth = min_depth;
        this.rmsMapQ = rmsMapQ;
        System.arraycopy(lk, 0, this.lk, 0, LIKELYHOOD_SIZE);
        byteSize = 9 + lk.length;
    }

    /**
     * the BinaryCodec constructor
     *
     * @param in the binary codec to get data from
     */
    SinglePointGenotype(BinaryCodec in) {
        offset = in.readInt();
        min_depth = in.readInt();
        rmsMapQ = in.readByte();
        in.readBytes(lk);
        byteSize = 9 + lk.length;
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
        out.writeBytes(lk);
    }

    public RECORD_TYPE getRecordType() {
        return RECORD_TYPE.SINGLE;
    }

}