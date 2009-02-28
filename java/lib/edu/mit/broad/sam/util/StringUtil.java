/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

public class StringUtil {
    /**
     *
     * @param separator String to interject between each string in strings arg
     * @param strings List of strings to be joined.
     * @return String that concatenates each item of strings arg, with separator btw each of them.
     */
    public static String join(final String separator, final String[] strings) {
        if (strings.length == 0) {
            return "";
        }
        final StringBuilder ret = new StringBuilder(strings[0]);
        for (int i = 1; i < strings.length; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }

    /**
     * Split the string into tokesn separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     *
     * @param aString  the string to split
     * @param tokens an array to hold the parsed tokens
     * @param delim  character that delimits tokens
     * @return the number of tokens parsed
     */
    public static int split(final String aString, final String[] tokens, final char delim) {

        final int maxTokens = tokens.length;
        int nTokens = 0;
        int start = 0;
        int end = aString.indexOf(delim);
        if(end < 0) {
            tokens[nTokens++] = aString;
            return nTokens;
        }
        while ((end > 0) && (nTokens < maxTokens))
        {
            tokens[nTokens++] = aString.substring(start, end);
            start = end + 1;
            end = aString.indexOf(delim, start);

        }
        // Add the trailing string,  if there is room and if it is not empty.
        if (nTokens < maxTokens)
        {
            final String trailingString = aString.substring(start);
            if (trailingString.length() > 0)
            {
                tokens[nTokens++] = trailingString;
            }
        }
        return nTokens;
    }

    ////////////////////////////////////////////////////////////////////
    // The following methods all convert btw bytes and Strings, without
    // using the Java character set mechanism.
    ////////////////////////////////////////////////////////////////////

    public static String bytesToString(final byte[] data) {
        if (data == null) {
            return null;
        }
        return bytesToString(data, 0, data.length);
    }

    @SuppressWarnings("deprecation")
    public static String bytesToString(final byte[] buffer, final int offset, final int length) {
/*
        The non-deprecated way, that requires allocating char[]
        final char[] charBuffer = new char[length];
        for (int i = 0; i < length; ++i) {
            charBuffer[i] = (char)buffer[i+offset];
        }
        return new String(charBuffer);
*/
        return new String(buffer, 0, offset, length);
    }

    @SuppressWarnings("deprecation")
    public static byte[] stringToBytes(final String s) {
/*
        The non-deprecated way, that requires allocating char[]
        final byte[] byteBuffer = new byte[s.length()];
        final char[] charBuffer = s.toCharArray();
        for (int i = 0; i < charBuffer.length; ++i) {
            byteBuffer[i] = (byte)(charBuffer[i] & 0xff);
        }
        return byteBuffer;
*/
        final byte[] byteBuffer = new byte[s.length()];
        s.getBytes(0, byteBuffer.length, byteBuffer, 0);
        return byteBuffer;
    }

    // This method might more appropriately live in BinaryCodec, but all the byte <=> char conversion
    // should be in the same place.
    public static String readNullTerminatedString(final BinaryCodec binaryCodec) {
        final StringBuilder ret = new StringBuilder();
        for (byte b = binaryCodec.readByte(); b != 0; b = binaryCodec.readByte()) {
            ret.append((char)(b & 0xff));
        }
        return ret.toString();
    }

    /**
     * Convert chars to bytes merely by casting
     * @param chars input chars
     * @param charOffset where to start converting from chars array
     * @param length how many chars to convert
     * @param bytes where to put the converted output
     * @param byteOffset where to start writing the converted output.
     */
    public static void charsToBytes(final char[] chars, final int charOffset, final int length,
                                    final byte[] bytes, final int byteOffset) {
        for (int i = 0; i < length; ++i) {
            bytes[byteOffset + i] = (byte)chars[charOffset + i];
        }
    }

}
