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

import java.io.File;
import java.util.List;
import java.util.Arrays;

/**
 * Parser for tab-delimited files
 *
 * @author Kathleen Tibbetts
 */
public class TabbedTextFileParser extends BasicTextFileParser {

    /**
     * Constructor
     *
     * @param file  The file to parse
     */
    public TabbedTextFileParser(boolean treatGroupedDelimitersAsOne, File... file) {
        super(treatGroupedDelimitersAsOne, file);
    }

    /**
     * Determines whether a given character is a delimiter
     *
     * @param b the character to evaluate
     * @return  true if <code>b</code> is a delimiter; otherwise false
     */
    protected boolean isDelimiter(byte b) {
        return b == '\t';
    }
}
