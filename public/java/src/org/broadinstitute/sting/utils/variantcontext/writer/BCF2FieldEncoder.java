/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext.writer;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Encoder;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Type;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCompoundHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 *
 *
 * @author Your Name
 * @since Date created
 */
@Invariant({
        "headerLine != null",
        "BCF2Type.INTEGERS.contains(dictionaryOffsetType)",
        "dictionaryOffset >= 0"
})
public abstract class BCF2FieldEncoder {
    final VCFCompoundHeaderLine headerLine;
    final BCF2Type fixedType;
    final int dictionaryOffset;
    final BCF2Type dictionaryOffsetType;

    // ----------------------------------------------------------------------
    //
    // Constructor
    //
    // ----------------------------------------------------------------------

    public BCF2FieldEncoder(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict, final BCF2Type fixedType) {
        this.headerLine = headerLine;
        this.fixedType = fixedType;

        final Integer offset = dict.get(getField());
        if ( offset == null ) throw new ReviewedStingException("Format error: could not find string " + getField() + " in header as required by BCF");
        this.dictionaryOffset = offset;
        dictionaryOffsetType = BCF2Utils.determineIntegerType(offset);
    }

    // ----------------------------------------------------------------------
    //
    // Basic accessors
    //
    // ----------------------------------------------------------------------

    public final String getField() { return headerLine.getID(); }

    /**
     * Write the field key (dictionary offset and type) into the BCF2Encoder stream
     *
     * @param encoder where we write our dictionary offset
     * @throws IOException
     */
    public final void writeFieldKey(final BCF2Encoder encoder) throws IOException {
        encoder.encodeTyped(dictionaryOffset, dictionaryOffsetType);
    }

    @Override
    public String toString() {
        return "BCF2FieldEncoder for " + getField() + " with count " + getCountType() + " encoded with " + getClass().getSimpleName();
    }

    // ----------------------------------------------------------------------
    //
    // methods to determine the number of encoded elements
    //
    // ----------------------------------------------------------------------

    protected final VCFHeaderLineCount getCountType() {
        return headerLine.getCountType();
    }

    @Ensures("result != (hasValueDeterminedNumElements() || hasContextDeterminedNumElements())")
    public boolean hasConstantNumElements() {
        return getCountType() == VCFHeaderLineCount.INTEGER;
    }

    @Ensures("result != (hasConstantNumElements() || hasContextDeterminedNumElements())")
    public boolean hasValueDeterminedNumElements() {
        return getCountType() == VCFHeaderLineCount.UNBOUNDED;
    }

    @Ensures("result != (hasValueDeterminedNumElements() || hasConstantNumElements())")
    public boolean hasContextDeterminedNumElements() {
        return ! hasConstantNumElements() && ! hasValueDeterminedNumElements();
    }

    @Requires("hasConstantNumElements()")
    @Ensures("result >= 0")
    public int numElements() {
        return headerLine.getCount();
    }

    @Requires("hasValueDeterminedNumElements()")
    @Ensures("result >= 0")
    public int numElements(final Object value) {
        return numElementsFromValue(value);
        //return value instanceof List ? ((List) value).size() : 1;
    }

    @Requires("hasContextDeterminedNumElements()")
    @Ensures("result >= 0")
    public int numElements(final VariantContext vc) {
        return headerLine.getCount(vc.getNAlleles() - 1);
    }

    @Ensures("result >= 0")
    public final int numElements(final VariantContext vc, final Object value) {
        if ( hasConstantNumElements() ) return numElements();
        else if ( hasContextDeterminedNumElements() ) return numElements(vc);
        else return numElements(value);
    }

    /**
     * Given a value, return the number of elements we will encode for it.
     *
     * Assumes the value is encoded as a List
     *
     * @param value
     * @return
     */
    @Requires("hasValueDeterminedNumElements()")
    @Ensures("result >= 0")
    protected int numElementsFromValue(final Object value) {
        if ( value == null ) return 0;
        else if ( value instanceof List ) return ((List) value).size();
        else return 1;
    }

    // ----------------------------------------------------------------------
    //
    // methods to determine the BCF2 type of the encoded values
    //
    // ----------------------------------------------------------------------

    @Ensures("result || isDynamicallyTyped()")
    public final boolean isStaticallyTyped() { return ! isDynamicallyTyped(); }

    @Ensures("result || isStaticallyTyped()")
    public final boolean isDynamicallyTyped() { return fixedType == null; }

    public final BCF2Type getType(final Object value) {
        return isDynamicallyTyped() ? getDynamicType(value) : getStaticType();
    }

    @Requires("isStaticallyTyped()")
    @Ensures("result != null")
    public final BCF2Type getStaticType() {
        return fixedType;
    }

    @Requires("isDynamicallyTyped()")
    @Ensures("result != null")
    public BCF2Type getDynamicType(final Object value) {
        throw new ReviewedStingException("BUG: cannot get dynamic type for statically typed BCF2 field");
    }

    // ----------------------------------------------------------------------
    //
    // methods to encode values, including the key abstract method
    //
    // ----------------------------------------------------------------------

    @Requires({"encoder != null", "isDynamicallyTyped() || type == getStaticType()"})
    public void encodeOneValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
        encodeValue(encoder, value, type, 0);
    }

    @Requires({"encoder != null", "isDynamicallyTyped() || type == getStaticType()", "minValues >= 0"})
    public abstract void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException;

    // ----------------------------------------------------------------------
    //
    // Subclass to encode Strings
    //
    // ----------------------------------------------------------------------

    public static class StringOrCharacter extends BCF2FieldEncoder {
        public StringOrCharacter(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.CHAR);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            final String s = javaStringToBCF2String(value);
            encoder.encodeString(s, Math.max(s.length(), minValues));
        }

        //
        // Regardless of what the header says, BCF2 strings and characters are always encoded
        // as arrays of CHAR type, which has a variable number of elements depending on the
        // exact string being encoded
        //
        @Override public boolean hasConstantNumElements()          { return false; }
        @Override public boolean hasContextDeterminedNumElements() { return false; }
        @Override public boolean hasValueDeterminedNumElements()   { return true; }
        @Override protected int numElementsFromValue(final Object value) {
            return value == null ? 0 : javaStringToBCF2String(value).length();
        }

        /**
         * Recode the incoming object to a String, compacting it into a
         * BCF2 string if the value is a list.
         *
         * @param value a String or List<String> to encode, or null
         * @return a non-null string to encode
         */
        @Ensures("result != null")
        private String javaStringToBCF2String(final Object value) {
            return value == null
                    ? ""
                    : (value instanceof List
                        ? BCF2Utils.collapseStringList((List<String>)value)
                        : (String)value);
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode FLAG
    //
    // ----------------------------------------------------------------------

    public static class Flag extends BCF2FieldEncoder {
        public Flag(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.INT8);
            if ( ! headerLine.isFixedCount() || headerLine.getCount() != 0 )
                throw new ReviewedStingException("Flag encoder only suppports atomic flags!");
        }

        @Override
        public int numElements() {
            return 1; // the header says 0 but we will write 1 value
        }

        @Override
        @Requires("minValues <= 1")
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            encoder.encodePrimitive(1, getStaticType());
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode FLOAT
    //
    // ----------------------------------------------------------------------

    public static class Float extends BCF2FieldEncoder {
        public Float(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.FLOAT);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            final List<Double> doubles = toList(Double.class, value);
            int count = 0;
            for ( final double d : doubles ) {
                encoder.encodeRawFloat(d);
                count++;
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode int[]
    //
    // ----------------------------------------------------------------------

    public static class IntArray extends BCF2FieldEncoder {
        public IntArray(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, null);
        }

        @Override
        protected int numElementsFromValue(final Object value) {
            return value == null ? 0 : ((int[])value).length;
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType((int[])value);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            if ( value != null ) {
                for ( final int i : (int[])value ) {
                    encoder.encodeRawInt(i, type);
                    count++;
                }
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode List<Integer>
    //
    // ----------------------------------------------------------------------

    public static class GenericInts extends BCF2FieldEncoder {
        public GenericInts(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType(toList(Integer.class, value));
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            for ( final int i : toList(Integer.class, value) ) {
                encoder.encodeRawInt(i, type);
                count++;
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }


    // ----------------------------------------------------------------------
    //
    // Helper methods
    //
    // ----------------------------------------------------------------------

    /**
     * Helper function that takes an object and returns a list representation
     * of it:
     *
     * o == null => []
     * o is a list => o
     * else => [o]
     *
     * @param o
     * @return
     */
    private final static <T> List<T> toList(final Class<T> c, final Object o) {
        if ( o == null ) return Collections.emptyList();
        else if ( o instanceof List ) return (List<T>)o;
        else return Collections.singletonList((T)o);
    }
}
