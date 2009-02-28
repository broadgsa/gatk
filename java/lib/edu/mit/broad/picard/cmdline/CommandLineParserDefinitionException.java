/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.cmdline;

public class CommandLineParserDefinitionException extends RuntimeException {
    public CommandLineParserDefinitionException() {
    }

    public CommandLineParserDefinitionException(String s) {
        super(s);
    }

    public CommandLineParserDefinitionException(String s, Throwable throwable) {
        super(s, throwable);
    }

    public CommandLineParserDefinitionException(Throwable throwable) {
        super(throwable);
    }
}
