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
 * Interface allows for implementations that read lines from a String, an ASCII file, or somewhere else.
 */
public interface LineReader {

    /**
     * Read a line and remove the line terminator
     */
    String readLine();

    /**
     * Read a line and optionally include the line terminator
     * @param includeTerminators
     * @return
     */
    String readLine(boolean includeTerminators);

    /**
     * @return 1-based number of line most recently read
     */
    int getLineNumber();
}
