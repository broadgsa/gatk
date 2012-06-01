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

package org.broadinstitute.sting.utils.codecs.bcf2;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

/**
 * Simple BCF2 encoder
 *
 * @author depristo
 * @since 5/12
 */
public class BCF2Encoder {
    // TODO -- increase default size?
    public static final int WRITE_BUFFER_INITIAL_SIZE = 16384;
    private ByteArrayOutputStream encodeStream = new ByteArrayOutputStream(WRITE_BUFFER_INITIAL_SIZE);

    // --------------------------------------------------------------------------------
    //
    // Functions to return the data being encoded here
    //
    // --------------------------------------------------------------------------------

    public int getRecordSizeInBytes() {
        return encodeStream.size();
    }

    public byte[] getRecordBytes() {
        byte[] bytes = encodeStream.toByteArray();
        encodeStream.reset();
        return bytes;
    }

    /**
     * Method for writing raw bytes to the encoder stream
     *
     * The purpuse this method exists is to allow lazy decoding of genotype data.  In that
     * situation the reader has loaded a block of bytes, and never decoded it, so we
     * are just writing it back out immediately as a raw stream of blocks.  Any
     * bad low-level formatting or changes to that byte[] will result in a malformed
     * BCF2 block.
     *
     * @param bytes a non-null byte array
     * @throws IOException
     */
    public void writeRawBytes(final byte[] bytes) throws IOException {
        assert bytes != null;
        encodeStream.write(bytes);
    }

    // --------------------------------------------------------------------------------
    //
    // Writing typed values (have type byte)
    //
    // --------------------------------------------------------------------------------

    public final void encodeTypedMissing(final BCF2Type type) throws IOException {
        encodeTyped(Collections.emptyList(), type);
    }

    // todo -- should be specialized for each object type for efficiency
    public final void encodeTyped(final Object v, final BCF2Type type) throws IOException {
        encodeTyped(Collections.singletonList(v), type);
    }

    public final void encodeTyped(List<? extends Object> v, final BCF2Type type) throws IOException {
        if ( type == BCF2Type.CHAR && v.size() != 0 ) {
            final String s = v.size() > 1 ? BCF2Utils.collapseStringList((List<String>)v) : (String)v.get(0);
            v = stringToBytes(s);
        }

        encodeType(v.size(), type);
        encodeRawValues(v, type);
    }

    // --------------------------------------------------------------------------------
    //
    // Writing raw values (don't have a type byte)
    //
    // --------------------------------------------------------------------------------

    public final <T extends Object> void encodeRawValues(final Collection<T> v, final BCF2Type type) throws IOException {
        for ( final T v1 : v ) {
            encodeRawValue(v1, type);
        }
    }

    public final <T extends Object> void encodeRawValue(final T value, final BCF2Type type) throws IOException {
        try {
            if ( value == type.getMissingJavaValue() )
                encodeRawMissingValue(type);
            else {
                switch (type) {
                    case INT8:
                    case INT16:
                    case INT32: encodePrimitive((Integer)value, type); break;
                    case FLOAT: encodeRawFloat((Double) value); break;
                    case CHAR:  encodeRawChar((Byte) value); break;
                    default:    throw new ReviewedStingException("Illegal type encountered " + type);
                }
            }
        } catch ( ClassCastException e ) {
            throw new ClassCastException("BUG: invalid type cast to " + type + " from " + value);
        }
    }

    public final void encodeRawMissingValue(final BCF2Type type) throws IOException {
        encodePrimitive(type.getMissingBytes(), type);
    }

    public final void encodeRawMissingValues(final int size, final BCF2Type type) throws IOException {
        if ( size <= 0 ) throw new ReviewedStingException("BUG: size <= 0");

        for ( int i = 0; i < size; i++ )
            encodeRawMissingValue(type);
    }

    // --------------------------------------------------------------------------------
    //
    // low-level encoders
    //
    // --------------------------------------------------------------------------------

    public final void encodeRawChar(final byte c) throws IOException {
        encodeStream.write(c);
    }

    public final void encodeRawFloat(final double value) throws IOException {
        encodePrimitive(Float.floatToIntBits((float)value), BCF2Type.FLOAT);
    }

    public final void encodeType(final int size, final BCF2Type type) throws IOException {
        if ( size < 0 ) throw new ReviewedStingException("BUG: size < 0");

        final byte typeByte = BCF2Utils.encodeTypeDescriptor(size, type);
        encodeStream.write(typeByte);
        if ( BCF2Utils.willOverflow(size) ) {
            // write in the overflow size
            encodeTyped(size, determineIntegerType(size));
        }
    }

    public final void encodeRawInt(final int value, final BCF2Type type) throws IOException {
        encodePrimitive(value, type, encodeStream);
    }

    public final void encodePrimitive(final int value, final BCF2Type type) throws IOException {
        encodePrimitive(value, type, encodeStream);
    }

    // --------------------------------------------------------------------------------
    //
    // utility functions
    //
    // --------------------------------------------------------------------------------

    public final BCF2Type determineIntegerType(final int[] values) {
        // literally a copy of the code below, but there's no general way to unify lists and arrays in java
        BCF2Type maxType = BCF2Type.INT8;
        for ( final int value : values ) {
            final BCF2Type type1 = determineIntegerType(value);
            switch ( type1 ) {
                case INT8: break;
                case INT16: maxType = BCF2Type.INT16; break;
                case INT32: return BCF2Type.INT32; // fast path for largest possible value
                default: throw new ReviewedStingException("Unexpected integer type " + type1 );
            }
        }
        return maxType;
    }

    public final BCF2Type determineIntegerType(final List<Integer> values) {
        BCF2Type maxType = BCF2Type.INT8;
        for ( final int value : values ) {
            final BCF2Type type1 = determineIntegerType(value);
            switch ( type1 ) {
                case INT8: break;
                case INT16: maxType = BCF2Type.INT16; break;
                case INT32: return BCF2Type.INT32; // fast path for largest possible value
                default: throw new ReviewedStingException("Unexpected integer type " + type1 );
            }
        }
        return maxType;
    }

    public final BCF2Type determineIntegerType(final int value) {
        for ( final BCF2Type potentialType : BCF2Utils.INTEGER_TYPES_BY_SIZE ) {
            if ( potentialType.withinRange(value) )
                return potentialType;
        }

        throw new ReviewedStingException("Integer cannot be encoded in allowable range of even INT32: " + value);
    }

    /**
     * Totally generic encoder that examines o, determines the best way to encode it, and encodes it
     *
     * This method is incredibly slow, but it's only used for UnitTests so it doesn't matter
     *
     * @param o
     * @return
     */
    protected final BCF2Type encode(final Object o) throws IOException {
        if ( o == null ) throw new ReviewedStingException("Generic encode cannot deal with null values");

        if ( o instanceof List ) {
            final BCF2Type type = determineBCFType(((List) o).get(0));
            encodeTyped((List) o, type);
            return type;
        } else {
            final BCF2Type type = determineBCFType(o);
            encodeTyped(o, type);
            return type;
        }
    }

    private final BCF2Type determineBCFType(final Object arg) {
        final Object toType = arg instanceof List ? ((List)arg).get(0) : arg;

        if ( toType instanceof Integer )
            return determineIntegerType((Integer)toType);
        else if ( toType instanceof String )
            return BCF2Type.CHAR;
        else if ( toType instanceof Double )
            return BCF2Type.FLOAT;
        else
            throw new ReviewedStingException("No native encoding for Object of type " + arg.getClass().getSimpleName());
    }

    public final static void encodePrimitive(final int value, final BCF2Type type, final OutputStream encodeStream) throws IOException {
        for ( int i = type.getSizeInBytes() - 1; i >= 0; i-- ) {
            final int shift = i * 8;
            int mask = 0xFF << shift;
            int byteValue = (mask & value) >> shift;
            encodeStream.write(byteValue);
        }
    }

    private final List<Byte> stringToBytes(final String v) throws IOException {
        if ( v == null || v.equals("") )
            return Collections.emptyList();
        else {
            // TODO -- this needs to be optimized away for efficiency
            final byte[] bytes = v.getBytes();
            final List<Byte> l = new ArrayList<Byte>(bytes.length);
            for ( int i = 0; i < bytes.length; i++) l.add(bytes[i]);
            return l;
        }
    }
}