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
abstract class GLFRecord {

    // fields common to all records
    protected REF_BASE refBase;
    protected long offset = 1;
    protected short minimumLikelihood = 0;
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

        private final short fieldValue;   // in kilograms

        /**
         * private constructor, used by the enum class to makes each enum value
         * @param value the short values specified in the enum listing
         */
        REF_BASE( short value ) {
            fieldValue = value;
        }

        /**
         * static method from returning a REF_BASE given the character representation
         *
         * @param value the character representation of a REF_BASE
         *
         * @return the corresponding REF_BASE
         * @throws IllegalArgumentException if the value passed can't be converted
         */
        public static REF_BASE toBase( char value ) {
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
    enum RECORD_TYPE {
        SINGLE((short) 1),
        VARIABLE((short) 2);

        private final short fieldValue;   // in kilograms

        RECORD_TYPE( short value ) {
            fieldValue = value;
        }

        public short getReadTypeValue() {
            return fieldValue;
        }
    }


    /**
     * Constructor, given the base a character reference base
     *
     * @param base              the reference base in the reference
     * @param offset            the offset from the beginning of the reference seq
     * @param minimumLikelihood it's minimum likelihood
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    public GLFRecord( char base, long offset, short minimumLikelihood, int readDepth, short rmsMapQ ) {
        REF_BASE newBase = REF_BASE.toBase(base);
        validateInput(newBase, offset, minimumLikelihood, readDepth, rmsMapQ);
    }

    /**
     * Constructor, given the base a REF_BASE
     *
     * @param base              the reference base in the reference
     * @param offset            the offset from the beginning of the reference seq
     * @param minimumLikelihood it's minimum likelihood
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    GLFRecord( REF_BASE base, long offset, short minimumLikelihood, int readDepth, short rmsMapQ ) {
        validateInput(base, offset, minimumLikelihood, readDepth, rmsMapQ);
    }

    /**
     * validate the input during construction, and store valid values
     *
     * @param base              the reference base in the reference, as a REF_BASE
     * @param offset            the offset from the beginning of the reference seq
     * @param minimumLikelihood it's minimum likelihood
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    private void validateInput( REF_BASE base, long offset, short minimumLikelihood, int readDepth, short rmsMapQ ) {
        this.refBase = base;
        if (offset > 4294967295L || offset < 0) {
            throw new IllegalArgumentException("Offset is out of bounds (0 to 0xffffffff) value passed = " + offset);
        }
        this.offset = offset;

        if (minimumLikelihood > 255 || minimumLikelihood < 0) {
            throw new IllegalArgumentException("minimumLikelihood is out of bounds (0 to 0xffffffff) value passed = " + minimumLikelihood);
        }
        this.minimumLikelihood = minimumLikelihood;

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
     void write( BinaryCodec out ) {
        out.writeUByte((short) ( this.getRecordType().getReadTypeValue() << 4 | ( refBase.getBaseHexValue() & 0x0f ) ));
        out.writeUInt(( (Long) offset ).intValue());
        out.writeUInt((new Long(readDepth).intValue()));
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
     * @param d a double value
     * @return a byte, capped at 
     */
    protected static short toCappedShort(double d) {
        return (d > 255.0) ? (byte)255 : (byte)Math.round(d);
    }

    /**
     * find the minimum value in a set of doubles
     * @param vals
     * @return
     */
    protected static double findMin(double vals[]) {
        if (vals.length < 1) throw new StingException("findMin: an array of size < 1 was passed in");

        double min = vals[0];
        for (double d: vals) {
            if (d < min) min = d;
        }
        return min;
    }

}

