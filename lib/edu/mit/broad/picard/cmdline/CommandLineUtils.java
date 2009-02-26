/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.cmdline;

import java.io.*;
import java.util.regex.Pattern;

public class CommandLineUtils {
    /** Regex for splitting on spaces. */
    public static final Pattern SPACE_SPLITTER = Pattern.compile(" ");

    // Regexes to split things apart on white space
    public static final Pattern TAB_SPLITTER = Pattern.compile("\\t");

    /** Checks that a file exists and is readable, and then returns a buffered reader for it. */
    public static BufferedReader getReader(File file) throws IOException {
        return new BufferedReader(new InputStreamReader(getInputStream(file)));
	}

    /** Checks that a file exists and is readable, and then returns a input stream  for it. */
    public static InputStream getInputStream(File file) throws IOException {
		if (!file.exists()) {
			throw new RuntimeException("Specified file does not exist: " + file);
		}

		if (!file.canRead()) {
			throw new RuntimeException("Specified file is not readable: " + file);
		}

		return new FileInputStream(file);
	}
}
