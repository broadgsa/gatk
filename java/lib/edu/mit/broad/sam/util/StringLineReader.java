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

/**
 * Implementation of LineReader that gets its input from a String.  No charset conversion
 * is necessary because the String is in unicode.  Handles CR, LF or CRLF line termination,
 * but if asked to return the line terminator, it always comes back as LF.
 */
public class StringLineReader implements LineReader {

    private final String theString;
    private int curPos = 0;
    private int lineNumber = 0;

    public StringLineReader(final String s) {
        // Simplify later processing by replacing crlf with just lf, and replacing solo cr with lf
        this.theString = s.replaceAll("\r\n", "\n").replaceAll("\r", "\n");
    }

    /**
     * Read a line and remove the line terminator
     */
    public String readLine() {
        return readLine(false);
    }

    /**
     * Read a line and optionally include the line terminator
     *
     * @param includeTerminators
     * @return
     */
    public String readLine(final boolean includeTerminators) {
        if (curPos == theString.length()) {
            return null;
        }
        final int nextLfIndex = theString.indexOf('\n', curPos);
        if (nextLfIndex == -1) {
            final int startPos = curPos;
            curPos = theString.length();
            ++lineNumber;
            return theString.substring(startPos);
        }
        final int startPos = curPos;
        final int endPos = nextLfIndex + (includeTerminators? 1: 0);
        curPos = nextLfIndex + 1;
        ++lineNumber;
        return theString.substring(startPos, endPos);
    }

    /**
     * @return 1-based number of line most recently read
     */
    public int getLineNumber() {
        return lineNumber;
    }
}
