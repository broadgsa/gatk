package org.broadinstitute.sting.utils.genotype.glf;

import net.sf.samtools.util.BinaryCodec;
import org.broadinstitute.sting.utils.StingException;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class RecordType
 *         <p/>
 *         The base record type for all GLF entries. Each record has a number of fields
 *         common to the record set.  This is also the source of the REF_BASE enumeration,
 *         which represents the accepted FASTA nucleotide symbols and their assocated GLF
 *         field values.
 */
public abstract class GLFRecord {
    public final static double LIKELIHOOD_SCALE_FACTOR = 10;


    // fields common to all records
    protected String contig;
    protected REF_BASE refBase;
    protected long position = 1;
    protected int readDepth = 0;
    protected short rmsMapQ = 0;

    // the size of this base structure
    protected final int baseSize = 10;

    /** the reference base enumeration, with their short (type) values for GLF */
    public enum REF_BASE {
        X((short) 0x00),
        A((short) 0x01),
        C((short) 0x02),
        M((short) 0x03),
        G((short) 0x04),
        R((short) 0x05),
        S((short) 0x06),
        V((short) 0x07),
        T((short) 0x08),
        W((short) 0x09),
        Y((short) 0x0A),
        H((short) 0x0B),
        K((short) 0x0C),
        D((short) 0x0D),
        B((short) 0x0E),
        N((short) 0x0F);

        private final short fieldValue;

        /**
         * private constructor, used by the enum class to makes each enum value
         *
         * @param value the short values specified in the enum listing
         */
        REF_BASE(short value) {
            fieldValue = value;
        }

        /**
         * return the character representation
         *
         * @return the char for the reference base
         */
        public char toChar() {
            return this.toString().charAt(0);
        }

        /**
         * static method from returning a REF_BASE given the character representation
         *
         * @param value the character representation of a REF_BASE
         *
         * @return the corresponding REF_BASE
         * @throws IllegalArgumentException if the value passed can't be converted
         */
        public static REF_BASE toBase(char value) {
            // for the case where they're passing in the enumated value
            if (value <= 0x0F && value >= 0) {
                return REF_BASE.values()[value];
            }
            String str = String.valueOf(value).toUpperCase();
            for (int x = 0; x < REF_BASE.values().length; x++) {
                if (REF_BASE.values()[x].toString().equals(str)) {
                    return REF_BASE.values()[x];
                }
            }
            throw new IllegalArgumentException("Counldn't find matching reference base for " + str);
        }

        /** @return the hex value of the given REF_BASE */
        public short getBaseHexValue() {
            return fieldValue;
        }
    }

    /** the record type enum, which enumerates the different records we can have in a GLF */
    public enum RECORD_TYPE {
        SINGLE((short) 1),
        VARIABLE((short) 2);

        private final short fieldValue;   // a place to store the type

        RECORD_TYPE(short value) {
            fieldValue = value;
        }

        public short getReadTypeValue() {
            return fieldValue;
        }
    }


    /**
     * Constructor, given the base a character reference base
     *
     * @param contig            the contig string
     * @param base              the reference base in the reference
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    public GLFRecord(String contig, char base, long position, int readDepth, short rmsMapQ) {
        REF_BASE newBase = REF_BASE.toBase(base);
        validateInput(contig, newBase, position, readDepth, rmsMapQ);
    }

    /**
     * Constructor, given the base a REF_BASE
     *
     * @param contig            the contig string
     * @param base              the reference base in the reference
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    GLFRecord(String contig, REF_BASE base, long position, int readDepth, short rmsMapQ) {
        validateInput(contig, base, position, readDepth, rmsMapQ);
    }

    /**
     * validate the input during construction, and store valid values
     *
     * @param chromosome        the reference contig, as a String
     * @param base              the reference base in the reference, as a REF_BASE
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    private void validateInput(String chromosome, REF_BASE base, long position, int readDepth, short rmsMapQ) {
        // add any validation to the contig string here
        this.contig = chromosome;

        this.refBase = base;

        if (position > 4294967295L || position < 0) {
            throw new IllegalArgumentException("Position is out of bounds (0 to 0xffffffff) value passed = " + position);
        }
        this.position = position;

//        if (minimumLikelihood > 255 || minimumLikelihood < 0) {
//            throw new IllegalArgumentException("minimumLikelihood is out of bounds (0 to 0xffffffff) value passed = " + minimumLikelihood);
//        }
//        this.minimumLikelihood = GLFRecord.toCappedShort(minimumLikelihood);

        if (readDepth > 16777215 || readDepth < 0) {
            throw new IllegalArgumentException("readDepth is out of bounds (0 to 0xffffff) value passed = " + readDepth);
        }
        this.readDepth = readDepth;

        if (rmsMapQ > 255 || rmsMapQ < 0) {
            throw new IllegalArgumentException("rmsMapQ is out of bounds (0 to 0xff) value passed = " + rmsMapQ);
        }
        this.rmsMapQ = rmsMapQ;
    }

    /**
     * write the this record to a binary codec output.
     *
     * @param out the binary codec to write to
     */
    void write(BinaryCodec out, GLFRecord lastRecord) {
        long offset = 0;
        if (lastRecord != null && lastRecord.getContig() == this.getContig())
            offset = this.position - lastRecord.getPosition();
        else
            offset = this.position - 1; // we start at one, we need to subtract that off
        short bite = ((short) (this.getRecordType().getReadTypeValue() << 4 | (refBase.getBaseHexValue() & 0x0f)));
        out.writeUByte((short) (this.getRecordType().getReadTypeValue() << 4 | (refBase.getBaseHexValue() & 0x0f)));
        out.writeUInt(((Long) (offset)).intValue()); // we have to subtract one, we're an offset
        long write = (long) ((long) (readDepth & 0xffffff) | (long) (this.getMinimumLikelihood() & 0xff) << 24);
        out.writeUInt(write);
        out.writeUByte((short) rmsMapQ);
    }

    /**
     * get the record type
     *
     * @return the record type enumeration
     */
    public abstract RECORD_TYPE getRecordType();

    /**
     * Return the size of this record in bytes.
     *
     * @return the size of this record type, in bytes
     */
    public int getByteSize() {
        return 10; // the record type field (1), offset (4), the min depth field (4), and the rms mapping (1)
    }

    /**
     * convert a double to a byte, capping it at the maximum value of 255
     *
     * @param d a double value
     *
     * @return a byte, capped at
     */
    protected static short toCappedShort(double d) {
        return (d > 255.0) ? (short) 255 : (short) Math.round(d);
    }

    /**
     * find the minimum value in a set of doubles
     *
     * @param vals
     *
     * @return
     */
    protected static double findMin(double vals[]) {
        if (vals.length < 1) throw new StingException("findMin: an array of size < 1 was passed in");

        double min = vals[0];
        for (double d : vals)
            if (d < min) min = d;

        return min;
    }

    public REF_BASE getRefBase() {
        return refBase;
    }

    public long getPosition() {
        return position;
    }

    public short getMinimumLikelihood() {
        return calculateMinLikelihood();
    }

    public int getReadDepth() {
        return readDepth;
    }

    public short getRmsMapQ() {
        return rmsMapQ;
    }

    public String getContig() {
        return this.contig;
    }

    /**
     * this method had to be abstracted so that the underlying records could set the minimum likelihood (ML) in the event
     * that the ML is above 255.  In this case the records need to scale the value appropriately, and warn the users.
     * @return a short of the minimum likelihood.
     */
    protected abstract short calculateMinLikelihood();
}

