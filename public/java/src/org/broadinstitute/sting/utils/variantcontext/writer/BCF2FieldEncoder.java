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
public abstract class BCF2FieldEncoder {
    final VCFCompoundHeaderLine headerLine;
    final BCF2Type fixedType;
    final int dictionaryOffset;
    final BCF2Type dictionaryOffsetType;

    public BCF2FieldEncoder(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict, final BCF2Type fixedType) {
        this.headerLine = headerLine;
        this.fixedType = fixedType;

        final Integer offset = dict.get(getField());
        if ( offset == null ) throw new ReviewedStingException("Format error: could not find string " + getField() + " in header as required by BCF");
        this.dictionaryOffset = offset;
        dictionaryOffsetType = BCF2Utils.determineIntegerType(offset);
    }

    public VCFHeaderLineCount getCountType() {
        return headerLine.getCountType();
    }

    public VCFCompoundHeaderLine getHeaderLine() {
        return headerLine;
    }

    public boolean hasFixedCount() { return getCountType() == VCFHeaderLineCount.INTEGER; }
    public boolean hasUnboundedCount() { return getCountType() == VCFHeaderLineCount.UNBOUNDED; }
    public boolean hasContextDeterminedCount() { return ! hasFixedCount() && ! hasUnboundedCount(); }

    @Requires("hasFixedCount()")
    public int getFixedCount() { return headerLine.getCount(); }
    public int getContextDeterminedCount(final VariantContext vc) {
        return headerLine.getCount(vc.getNAlleles() - 1);
    }
    public int getBCFFieldCount(final VariantContext vc, final Object value) {
        if ( hasFixedCount() )
            return getFixedCount();
        else if ( hasUnboundedCount() )
            return value instanceof List ? ((List) value).size() : 1;
        else
            return getContextDeterminedCount(vc);
    }

    public String getField() { return headerLine.getID(); }

    public int getDictionaryOffset() { return dictionaryOffset; }
    public BCF2Type getDictionaryOffsetType() { return dictionaryOffsetType; }

    public boolean isFixedTyped() { return ! isDynamicallyTyped(); }
    public boolean isDynamicallyTyped() { return fixedType == null; }
    public BCF2Type getType(final Object value) { return isDynamicallyTyped() ? getDynamicType(value) : getFixedType(); }
    public BCF2Type getFixedType() {
        if ( fixedType != null )
            return fixedType;
        else
            throw new ReviewedStingException("Not a fixed type encoder: " + getField());
    }
    public BCF2Type getDynamicType(final Object value) { throw new ReviewedStingException("Function getDynamicType() not implemented"); }

    @Override
    public String toString() {
        return "BCF2FieldEncoder for " + getField() + " with count " + getCountType() + " encoded with " + getClass().getSimpleName();
    }

    public abstract void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException;


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

    public static class StringOrCharacter extends BCF2FieldEncoder {
        public StringOrCharacter(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, BCF2Type.CHAR);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            if ( value != null ) {
                final String s = encodeString(value);
                encoder.encodeString(s, s.length());
            }
        }

        @Override
        public int getBCFFieldCount(final VariantContext vc, final Object value) {
            return value == null ? 0 : encodeString(value).length();
        }

        private String encodeString(final Object value) {
            return value instanceof List ? BCF2Utils.collapseStringList((List<String>)value) : (String)value;
        }
    }

    public static class Flag extends BCF2FieldEncoder {
        public Flag(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, BCF2Type.INT8);
            if ( getHeaderLine().getCount() != 0 )
                throw new ReviewedStingException("Flag encoder only suppports atomic flags!");
        }

        @Override
        public int getFixedCount() {
            return 1; // the header says 0 but we will write 1 value
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            encoder.encodePrimitive(1, getFixedType());
        }
    }

    public static class Float extends BCF2FieldEncoder {
        public Float(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, BCF2Type.FLOAT);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            final List<Double> doubles = toList(Double.class, value);
            for ( final double d : doubles )
                encoder.encodeRawFloat(d);
        }
    }

    public static class IntArray extends BCF2FieldEncoder {
        public IntArray(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType((int[])value);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            for ( final int i : (int[])value )
                encoder.encodeRawInt(i, type);
        }
    }

    public static class IntList extends BCF2FieldEncoder {
        public IntList(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType(toList(Integer.class, value));
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            for ( final int i : toList(Integer.class, value) )
                encoder.encodeRawInt(i, type);
        }
    }

    public static class AtomicInt extends BCF2FieldEncoder {
        public AtomicInt(final VCFCompoundHeaderLine headerLine, final BCF2Encoder encoder, final Map<String, Integer> dict ) {
            super(headerLine, encoder, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType((Integer)value);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type) throws IOException {
            encoder.encodeRawInt((Integer)value, type);
        }
    }
}
