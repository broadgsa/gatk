/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.RuntimeEOFException;

import java.util.Map;
import java.util.Collection;
import java.util.ArrayList;

/**
 * Parse & produce tag section of alignment record in BAM file.
 */
class BinaryTagCodec {
    // Size of the fixed part of the binary representation of a tag,
    // i.e. the number of bytes occupied by the tag name and tag type fields.
    private static final int FIXED_TAG_SIZE = 3;

    private static final long MAX_INT = Integer.MAX_VALUE;
    private static final long MAX_UINT = (MAX_INT + 1) * 2;
    private static final long MAX_SHORT = Short.MAX_VALUE;
    private static final long MAX_USHORT = (MAX_SHORT + 1) * 2;
    private static final long MAX_BYTE = Byte.MAX_VALUE;
    private static final long MAX_UBYTE = (MAX_BYTE + 1) * 2;

    final BinaryCodec binaryCodec;

    BinaryTagCodec(final BinaryCodec binaryCodec) {
        this.binaryCodec = binaryCodec;
    }

    private static int getBinaryValueSize(final Object attributeValue) {
        switch (getTagValueType(attributeValue)) {
            case 'Z':
                return ((String)attributeValue).length() + 1;
            case 'A':
                return 1;
            case 'I':
            case 'i':
                return 4;
            case 's':
            case 'S':
                return 2;
            case 'c':
            case 'C':
                return 1;
            case 'f':
                return 4;
            case 'H':
                final byte[] byteArray = (byte[])attributeValue;
                return byteArray.length * 2 + 1;
            default:
                throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                        attributeValue.getClass().getName());
        }
    }

    static int getTagSize(final Object value) {
        return FIXED_TAG_SIZE + getBinaryValueSize(value);
    }

    static char getTagValueType(final Object value) {
        if (value.getClass().equals(String.class)) {
            return 'Z';
        } else if (value.getClass().equals(Character.class)) {
            return 'A';
        } else if (value.getClass().equals(Integer.class)) {
            return getIntegerType((Integer)value);
        } else if (value.getClass().equals(Long.class)) {
            return getIntegerType((Long)value);
        } else if (value.getClass().equals(Float.class)) {
            return 'f';
        } else if (value.getClass().isArray() && value.getClass().getComponentType().equals(Byte.class)) {
            return 'H';
        } else {
            throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                    value.getClass().getName());
        }
    }

    static private char getIntegerType(final long val) {
        if (val > MAX_UINT) {
            throw new IllegalArgumentException("Integer attribute value too large to be encoded in BAM");
        }
        if (val > MAX_INT) {
            return 'I';
        }
        if (val > MAX_USHORT) {
            return 'i';
        }
        if (val > MAX_SHORT) {
            return 'S';
        }
        if (val > MAX_UBYTE) {
            return 's';
        }
        if (val > MAX_BYTE) {
            return 'C';
        }
        if (val >= Byte.MIN_VALUE) {
            return 'c';
        }
        if (val >= Short.MIN_VALUE) {
            return 's';
        }
        if (val >= Integer.MIN_VALUE) {
            return 'i';
        }
        throw new IllegalArgumentException("Integer attribute value too negative to be encoded in BAM");
    }

    void writeTag(final String key, final Object value) {
        assert(key.length() == 2);
        binaryCodec.writeString(key, false, false);
        final char tagValueType = getTagValueType(value);
        binaryCodec.writeByte(tagValueType);

        switch (tagValueType) {
            case 'Z':
                binaryCodec.writeString((String)value, false, true);
                break;
            case 'A':
                binaryCodec.writeByte(((Character)value));
                break;
            case 'I':
                binaryCodec.writeUInt((Long)value);
                break;
            case 'i':
                binaryCodec.writeInt((Integer)value);
                break;
            case 's':
                binaryCodec.writeShort(((Integer)value).shortValue());
                break;
            case 'S':
                binaryCodec.writeUShort((Integer)value);
                break;
            case 'c':
                binaryCodec.writeByte((Integer)value);
                break;
            case 'C':
                binaryCodec.writeUByte(((Integer)value).shortValue());
                break;
            case 'f':
                binaryCodec.writeFloat((Float)value);
                break;
            case 'H':
                final byte[] byteArray = (byte[])value;
                binaryCodec.writeString(SAMUtils.bytesToHexString(byteArray), false, true);
                break;
            default:
                throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                        value.getClass().getName());
        }
    }

    /**
     * Reads tags from the binaryCodec passed in the ctor
     * @param tagCollection tags are stored in this Map
     */
    void readTags(final Map<String, Object> tagCollection) {
        while (true) {
            final String key;
            try {
                // Only way to know at end is when out of input
                key = binaryCodec.readString(2);
            } catch (RuntimeEOFException e) {
                break;
            }
            final byte tagType = binaryCodec.readByte();
            final Object value = readValue(tagType);
            tagCollection.put(key, value);
        }
    }

    private Object readValue(final byte tagType) {
        switch (tagType) {
            case 'Z':
                return binaryCodec.readNullTerminatedString();
            case 'A':
                return (char)binaryCodec.readByte();
            case 'I':
                return binaryCodec.readUInt();
            case 'i':
                return binaryCodec.readInt();
            case 's':
                return (int)binaryCodec.readShort();
            case 'S':
                return binaryCodec.readUShort();
            case 'c':
                return (int)binaryCodec.readByte();
            case 'C':
                return (int)binaryCodec.readUByte();
            case 'f':
                return binaryCodec.readFloat();
            case 'H':
                final String hexRep = binaryCodec.readNullTerminatedString();
                return SAMUtils.hexStringToBytes(hexRep);
            default:
                throw new SAMFormatException("Unrecognized tag type: " + (char)tagType);
        }
    }

}
