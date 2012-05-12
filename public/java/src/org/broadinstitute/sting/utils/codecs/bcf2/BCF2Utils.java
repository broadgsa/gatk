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

import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFIDHeaderLine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Common utilities for working with BCF2 files
 *
 * Includes convenience methods for encoding, decoding BCF2 type descriptors (size + type)
 *
 * @author depristo
 * @since 5/12
 */
public class BCF2Utils {
    public static final byte[] MAGIC_HEADER_LINE = "BCF\2".getBytes();

    public static final int OVERFLOW_ELEMENT_MARKER = 15;
    public static final int MAX_INLINE_ELEMENTS = 14;

    // Note that these values are prefixed by FFFFFF for convenience
    public static final int INT8_MISSING_VALUE  = 0xFFFFFF80;
    public static final int INT16_MISSING_VALUE = 0xFFFF8000;
    public static final int INT32_MISSING_VALUE = 0x80000000;
    public static final int FLOAT_MISSING_VALUE = 0x7F800001;

    public final static BCF2Type[] INTEGER_TYPES_BY_SIZE = new BCF2Type[]{BCF2Type.INT8, BCF2Type.INT16, BCF2Type.INT32};

    private BCF2Utils() {}

    /**
     * Create a strings dictionary from the VCF header
     *
     * The dictionary is an ordered list of common VCF identifers (FILTER, INFO, and FORMAT)
     * fields.
     *
     * @param header the VCFHeader from which to build the dictionary
     * @return a non-null dictionary of elements, may be empty
     */
    public final static ArrayList<String> makeDictionary(final VCFHeader header) {
        final ArrayList<String> dict = new ArrayList<String>();

        // set up the strings dictionary
        dict.add(VCFConstants.PASSES_FILTERS_v4); // special case the special PASS field
        for ( VCFHeaderLine line : header.getMetaData() ) {
            if ( line instanceof VCFIDHeaderLine) {
                VCFIDHeaderLine idLine = (VCFIDHeaderLine)line;
                dict.add(idLine.getID());
            }
        }

        return dict;
    }

    public final static byte encodeTypeDescriptor(final int nElements, final BCF2Type type ) {
        int encodeSize = Math.min(nElements, OVERFLOW_ELEMENT_MARKER);
        byte typeByte = (byte)((0x0F & encodeSize) << 4 | (type.getID() & 0x0F));
        return typeByte;
    }

    public final static int decodeSize(final byte typeDescriptor) {
        return (0xF0 & typeDescriptor) >> 4;
    }

    public final static int decodeTypeID(final byte typeDescriptor) {
        return typeDescriptor & 0x0F;
    }

    public final static BCF2Type decodeType(final byte typeDescriptor) {
        return BCF2Type.values()[decodeTypeID(typeDescriptor)];
    }

    public final static boolean sizeIsOverflow(final byte typeDescriptor) {
        return decodeSize(typeDescriptor) == OVERFLOW_ELEMENT_MARKER;
    }

    public final static boolean willOverflow(final long nElements) {
        return nElements > MAX_INLINE_ELEMENTS;
    }

    public final static boolean startsWithBCF2Magic(final InputStream stream) throws IOException {
        final byte[] magicBytes = new byte[BCF2Utils.MAGIC_HEADER_LINE.length];
        stream.read(magicBytes);
        return Arrays.equals(magicBytes, BCF2Utils.MAGIC_HEADER_LINE);
    }

    public final static byte readByte(final InputStream stream) {
        try {
            return (byte)(stream.read() & 0xFF);
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }

    public final static int readInt(int bytesForEachInt, final InputStream stream) {
        switch ( bytesForEachInt ) {
            case 1: {
                return (byte)(readByte(stream));
            } case 2: {
                final int b1 = readByte(stream) & 0xFF;
                final int b2 = readByte(stream) & 0xFF;
                return (short)((b1 << 8) | b2);
            } case 4: {
                final int b1 = readByte(stream) & 0xFF;
                final int b2 = readByte(stream) & 0xFF;
                final int b3 = readByte(stream) & 0xFF;
                final int b4 = readByte(stream) & 0xFF;
                return (int)(b1 << 24 | b2 << 16 | b3 << 8 | b4);
            } default: throw new ReviewedStingException("Unexpected size during decoding");
        }
    }
}
