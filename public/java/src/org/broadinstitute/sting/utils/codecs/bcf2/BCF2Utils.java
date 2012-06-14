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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFIDHeaderLine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;

/**
 * Common utilities for working with BCF2 files
 *
 * Includes convenience methods for encoding, decoding BCF2 type descriptors (size + type)
 *
 * @author depristo
 * @since 5/12
 */
public final class BCF2Utils {
    public static final byte[] MAGIC_HEADER_LINE = "BCF\2".getBytes();

    public static final int MAX_ALLELES_IN_GENOTYPES = 127;

    public static final int OVERFLOW_ELEMENT_MARKER = 15;
    public static final int MAX_INLINE_ELEMENTS = 14;

    public final static BCF2Type[] INTEGER_TYPES_BY_SIZE = new BCF2Type[]{BCF2Type.INT8, BCF2Type.INT16, BCF2Type.INT32};
    public final static BCF2Type[] ID_TO_ENUM;

    static {
        int maxID = -1;
        for ( BCF2Type v : BCF2Type.values() ) maxID = Math.max(v.getID(), maxID);
        ID_TO_ENUM = new BCF2Type[maxID+1];
        for ( BCF2Type v : BCF2Type.values() ) ID_TO_ENUM[v.getID()] = v;
    }

    private BCF2Utils() {}

    /**
     * Create a strings dictionary from the VCF header
     *
     * The dictionary is an ordered list of common VCF identifers (FILTER, INFO, and FORMAT)
     * fields.
     *
     * Note that its critical that the list be dedupped and sorted in a consistent manner each time,
     * as the BCF2 offsets are encoded relative to this dictionary, and if it isn't determined exactly
     * the same way as in the header each time it's very bad
     *
     * @param header the VCFHeader from which to build the dictionary
     * @return a non-null dictionary of elements, may be empty
     */
    @Requires("header != null")
    @Ensures({"result != null", "new HashSet(result).size() == result.size()"})
    public final static ArrayList<String> makeDictionary(final VCFHeader header) {
        final Set<String> dict = new TreeSet<String>();

        // set up the strings dictionary
        dict.add(VCFConstants.PASSES_FILTERS_v4); // special case the special PASS field
        for ( VCFHeaderLine line : header.getMetaData() ) {
            if ( line instanceof VCFIDHeaderLine) {
                VCFIDHeaderLine idLine = (VCFIDHeaderLine)line;
                dict.add(idLine.getID());
            }
        }

        return new ArrayList<String>(dict);
    }

    @Requires({"nElements >= 0", "type != null"})
    public final static byte encodeTypeDescriptor(final int nElements, final BCF2Type type ) {
        int encodeSize = Math.min(nElements, OVERFLOW_ELEMENT_MARKER);
        byte typeByte = (byte)((0x0F & encodeSize) << 4 | (type.getID() & 0x0F));
        return typeByte;
    }

    @Ensures("result >= 0")
    public final static int decodeSize(final byte typeDescriptor) {
        return (0xF0 & typeDescriptor) >> 4;
    }

    @Ensures("result >= 0")
    public final static int decodeTypeID(final byte typeDescriptor) {
        return typeDescriptor & 0x0F;
    }

    @Ensures("result != null")
    public final static BCF2Type decodeType(final byte typeDescriptor) {
        return ID_TO_ENUM[decodeTypeID(typeDescriptor)];
    }

    public final static boolean sizeIsOverflow(final byte typeDescriptor) {
        return decodeSize(typeDescriptor) == OVERFLOW_ELEMENT_MARKER;
    }

    @Requires("nElements >= 0")
    public final static boolean willOverflow(final long nElements) {
        return nElements > MAX_INLINE_ELEMENTS;
    }

    public final static boolean startsWithBCF2Magic(final InputStream stream) throws IOException {
        final byte[] magicBytes = new byte[BCF2Utils.MAGIC_HEADER_LINE.length];
        stream.read(magicBytes);
        return Arrays.equals(magicBytes, BCF2Utils.MAGIC_HEADER_LINE);
    }

    public final static byte readByte(final InputStream stream) {
        // TODO -- shouldn't be capturing error here
        try {
            return (byte)(stream.read() & 0xFF);
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }

    @Requires({"stream != null", "bytesForEachInt > 0"})
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

    /**
     * Collapse multiple strings into a comma separated list
     *
     * ["s1", "s2", "s3"] => ",s1,s2,s3"
     *
     * @param strings size > 1 list of strings
     * @return
     */
    @Requires({"strings != null", "strings.size() > 1"})
    @Ensures("result != null")
    public static final String collapseStringList(final List<String> strings) {
        final StringBuilder b = new StringBuilder();
        for ( final String s : strings ) {
            assert s.indexOf(",") == -1; // no commas in individual strings
            b.append(",").append(s);
        }
        return b.toString();
    }

    /**
     * Inverse operation of collapseStringList.
     *
     * ",s1,s2,s3" => ["s1", "s2", "s3"]
     *
     *
     * @param collapsed
     * @return
     */
    @Requires({"collapsed != null", "isCollapsedString(collapsed)"})
    @Ensures("result != null")
    public static final List<String> exploreStringList(final String collapsed) {
        assert isCollapsedString(collapsed);
        final String[] exploded = collapsed.substring(1).split(",");
        return Arrays.asList(exploded);
    }

    @Requires("s != null")
    public static final boolean isCollapsedString(final String s) {
        return s.charAt(0) == ',';
    }

    /**
     * Returns a good name for a shadow BCF file for vcfFile.
     *
     * foo.vcf => foo.bcf
     * foo.xxx => foo.xxx.bcf
     *
     * @param vcfFile
     * @return
     */
    @Requires("vcfFile != null")
    @Ensures("result != null")
    public static final File shadowBCF(final File vcfFile) {
        final String path = vcfFile.getAbsolutePath();
        if ( path.contains(".vcf") )
            return new File(path.replace(".vcf", ".bcf"));
        else
            return new File( path + ".bcf" );
    }

    @Ensures("BCF2Type.INTEGERS.contains(result)")
    public final static BCF2Type determineIntegerType(final int value) {
        for ( final BCF2Type potentialType : INTEGER_TYPES_BY_SIZE) {
            if ( potentialType.withinRange(value) )
                return potentialType;
        }

        throw new ReviewedStingException("Integer cannot be encoded in allowable range of even INT32: " + value);
    }

    @Ensures("BCF2Type.INTEGERS.contains(result)")
    public final static BCF2Type determineIntegerType(final int[] values) {
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

    /**
     * Returns the maximum BCF2 integer size of t1 and t2
     *
     * For example, if t1 == INT8 and t2 == INT16 returns INT16
     *
     * @param t1
     * @param t2
     * @return
     */
    @Requires({"BCF2Type.INTEGERS.contains(t1)","BCF2Type.INTEGERS.contains(t2)"})
    @Ensures("BCF2Type.INTEGERS.contains(result)")
    public final static BCF2Type maxIntegerType(final BCF2Type t1, final BCF2Type t2) {
        switch ( t1 ) {
            case INT8: return t2;
            case INT16: return t2 == BCF2Type.INT32 ? t2 : t1;
            case INT32: return t1;
            default: throw new ReviewedStingException("BUG: unexpected BCF2Type " + t1);
        }
    }

    @Ensures("BCF2Type.INTEGERS.contains(result)")
    public final static BCF2Type determineIntegerType(final List<Integer> values) {
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
    public final static List<Object> toList(final Object o) {
        if ( o == null ) return Collections.emptyList();
        else if ( o instanceof List ) return (List<Object>)o;
        else return Collections.singletonList(o);
    }

    public final static void encodeRawBytes(final int value, final BCF2Type type, final OutputStream encodeStream) throws IOException {
        switch ( type.getSizeInBytes() ) {
            case 1:
                encodeStream.write(0xFF & value);
                break;
            case 2:
                encodeStream.write((0xFF00 & value) >> 8);
                encodeStream.write(0xFF & value);
                break;
            case 4:
                encodeStream.write((0xFF000000 & value) >> 24);
                encodeStream.write((0x00FF0000 & value) >> 16);
                encodeStream.write((0x0000FF00 & value) >> 8);
                encodeStream.write((0x000000FF & value));
                break;
            default:
                throw new ReviewedStingException("BUG: unexpected type size " + type);
        }
// general case for reference
//        for ( int i = type.getSizeInBytes() - 1; i >= 0; i-- ) {
//            final int shift = i * 8;
//            int mask = 0xFF << shift;
//            int byteValue = (mask & value) >> shift;
//            encodeStream.write(byteValue);
//        }
    }
}
