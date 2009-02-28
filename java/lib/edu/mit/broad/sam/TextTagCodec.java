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

import edu.mit.broad.sam.util.StringUtil;

import java.util.Map;

class TextTagCodec {
    private static final int NUM_TAG_FIELDS = 3;

    /**
     * This is really a local variable of decode(), but allocated here to reduce allocations.
     */
    private final String[] fields = new String[NUM_TAG_FIELDS];

    String encode(final String key, Object value) {
        final StringBuilder sb = new StringBuilder(key);
        sb.append(':');
        char tagType = BinaryTagCodec.getTagValueType(value);
        switch (tagType) {
            case 'c':
            case 'C':
            case 's':
            case 'S':
            case 'I':
                tagType = 'i';
        }
        if (tagType == 'H') {
            value = SAMUtils.bytesToHexString((byte[])value);
        }
        sb.append(tagType);
        sb.append(':');
        sb.append(value.toString());
        return sb.toString();
    }

    Map.Entry<String, Object> decode(final String tag) {
        final int numFields = StringUtil.split(tag, fields, ':');
        if (numFields != TextTagCodec.NUM_TAG_FIELDS) {
            throw new SAMFormatException("Not enough fields in tag '" + tag + "'");
        }
        final String key = fields[0];
        final String type = fields[1];
        final String stringVal = fields[2];
        final Object val;
        if (type.equals("Z")) {
            val = stringVal;
        } else if (type.equals("A")) {
            if (stringVal.length() != 1) {
                throw new SAMFormatException("Tag of type A should have a single-character value");
            }
            val = stringVal.charAt(0);
        } else if (type.equals("i")) {
            try {
                val = new Integer(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type i should have signed decimal value");
            }
        } else if (type.equals("f")) {
            try {
                val = new Float(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type f should have single-precision floating point value");
            }
        } else if (type.equals("H")) {
            try {
                val = SAMUtils.hexStringToBytes(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type H should have valid hex string with even number of digits");
            }
        } else {
            throw new SAMFormatException("Unrecognized tag type: " + type);
        }
        return new Map.Entry<String, Object>() {
            public String getKey() {
                return key;
            }

            public Object getValue() {
                return val;
            }

            public Object setValue(final Object o) {
                throw new UnsupportedOperationException();
            }
        };
    }
}
