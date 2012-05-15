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

    // --------------------------------------------------------------------------------
    //
    // Super-high level interface
    //
    // --------------------------------------------------------------------------------

    /**
     * Totally generic encoder that examines o, determines the best way to encode it, and encodes it
     * @param o
     * @return
     */
    public final BCF2Type encode(final Object o) throws IOException {
        if ( o == null ) throw new ReviewedStingException("Generic encode cannot deal with null values");

        if ( o instanceof String ) {
            return encodeString((String)o);
        } else if ( o instanceof List ) {
            final BCF2Type type = determinePrimitiveType(((List) o).get(0));
            encodeTypedVector((List) o, type);
            return type;
        } else {
            final BCF2Type type = determinePrimitiveType(o);
            encodeTypedSingleton(o, type);
            return type;
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Writing typed values (have type byte)
    //
    // --------------------------------------------------------------------------------

    public final void encodeTypedMissing(final BCF2Type type) throws IOException {
        encodeTypedVector(Collections.emptyList(), type);
    }

    // todo -- should be specialized for each object type for efficiency
    public final void encodeTypedSingleton(final Object v, final BCF2Type type) throws IOException {
        encodeTypedVector(Collections.singleton(v), type);
    }

    public final BCF2Type encodeString(final String v) throws IOException {
        // TODO -- this needs to be optimized
        final byte[] bytes = v.getBytes();
        final List<Byte> l = new ArrayList<Byte>(bytes.length);
        for ( int i = 0; i < bytes.length; i++) l.add(bytes[i]);
        encodeTypedVector(l, BCF2Type.CHAR);
        return BCF2Type.CHAR;
    }

    public final <T extends Object> void encodeTypedVector(final Collection<T> v, final BCF2Type type) throws IOException {
        encodeType(v.size(), type);
        encodeRawValues(v, type);
    }

    public final BCF2Type encodeTypedIntOfBestSize(final int value) throws IOException {
        final BCF2Type type = determineIntegerType(value);
        encodeTypedSingleton(value, type);
        return type;
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
        if ( value == type.getMissingJavaValue() )
            encodeRawMissingValue(type);
        else {
            switch (type) {
                case INT8:
                case INT16:
                case INT32: encodePrimitive((Integer)value, type); break;
                case FLOAT: encodeRawFloat((Double) value, type); break;
                case CHAR:  encodeRawChar((Byte) value); break;
                default:    throw new ReviewedStingException("Illegal type encountered " + type);
            }
        }
    }

    public final void encodeRawMissingValue(final BCF2Type type) throws IOException {
        encodePrimitive(type.getMissingBytes(), type);
    }

    public final void encodeRawMissingValues(final int size, final BCF2Type type) throws IOException {
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

    public final void encodeRawFloat(final double value, final BCF2Type type) throws IOException {
        encodePrimitive(Float.floatToIntBits((float)value), type);
    }

    public final void encodeType(final int size, final BCF2Type type) throws IOException {
        final byte typeByte = BCF2Utils.encodeTypeDescriptor(size, type);
        encodeStream.write(typeByte);
        if ( BCF2Utils.willOverflow(size) )
            encodeTypedIntOfBestSize(size);
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

    private final BCF2Type determinePrimitiveType(final Object v) {
        if ( v instanceof Integer )
            return determineIntegerType((Integer)v);
        else if ( v instanceof Double )
            return BCF2Type.FLOAT;
        else
            throw new ReviewedStingException("No native encoding for Object of type " + v.getClass().getSimpleName());
    }

    public final static void encodePrimitive(final int value, final BCF2Type type, final OutputStream encodeStream) throws IOException {
        for ( int i = type.getSizeInBytes() - 1; i >= 0; i-- ) {
            final int shift = i * 8;
            int mask = 0xFF << shift;
            int byteValue = (mask & value) >> shift;
            encodeStream.write(byteValue);
        }
    }
}
