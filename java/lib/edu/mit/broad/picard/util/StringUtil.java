/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

/**
 * Utilities that are useful when dealing with Strings.
 *
 * @author Tim Fennell
 */
public class StringUtil {
    /**
     * Return input string with newlines inserted to ensure that all lines
     * have length <= maxLineLength.  if a word is too long, it is simply broken
     * at maxLineLength.  Does not handle tabs intelligently (due to implementer laziness).
     */
    public static String wordWrap(String s, int maxLineLength) {
        String[] lines = s.split("\n");
        StringBuilder sb = new StringBuilder();
        for (String line: lines) {
            if (sb.length() > 0) {
                sb.append("\n");
            }
            sb.append(wordWrapSingleLine(line, maxLineLength));
        }
        if (s.endsWith("\n")) {
            sb.append("\n");
        }
        return sb.toString();
    }

    public static String wordWrapSingleLine(String s, int maxLineLength) {
        if (s.length() <= maxLineLength) {
            return s;
        }
        StringBuilder sb = new StringBuilder();
        int startCopyFrom = 0;
        while (startCopyFrom < s.length()) {
            int lastSpaceIndex = startCopyFrom;
            int i;
            // Find break point (if it exists)
            for (i = startCopyFrom; i < s.length() && i - startCopyFrom < maxLineLength; ++i) {
                if (Character.isWhitespace(s.charAt(i))) {
                    lastSpaceIndex = i;
                }
            }
            if (i - startCopyFrom < maxLineLength) {
                lastSpaceIndex = i;
            }
            // Include any trailing whitespace
            for (; lastSpaceIndex < s.length() && Character.isWhitespace(s.charAt(lastSpaceIndex)); ++lastSpaceIndex) {}
            if (sb.length() > 0) {
                sb.append("\n");
            }
            // Handle situation in which there is no word break.  Just break the word in the middle.
            if (lastSpaceIndex == startCopyFrom) {
                lastSpaceIndex = i;
            }
            sb.append(s.substring(startCopyFrom, lastSpaceIndex));
            startCopyFrom = lastSpaceIndex;
        }
        return sb.toString();
    }

    /**
     *
     * @param separator String to interject between each string in strings arg
     * @param strings List of strings to be joined.
     * @return String that concatenates each item of strings arg, with separator btw each of them.
     */
    public static String join(String separator, String... strings) {
        if (strings.length == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[0]);
        for (int i = 1; i < strings.length; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }

    /**
     * Checks that a String doesn't contain one or more characters of interest.
     *
     * @param s the String to check
     * @param chars the characters to check for
     * @return String the input String for convenience
     * @throws IllegalArgumentException if the String contains one or more of the characters
     */
    public static String assertCharactersNotInString(final String s, final char... chars) {
        for (char ch : s.toCharArray()) {
            for (int i=0; i<chars.length; ++i) {
                if (ch == chars[i]) {
                    throw new IllegalArgumentException("Supplied String contains illegal character '" + chars[i] + "'.");
                }
            }
        }

        return s;
    }
}
