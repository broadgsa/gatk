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

public class RuntimeEOFException extends RuntimeException {
    public RuntimeEOFException() {
    }

    public RuntimeEOFException(final String s) {
        super(s);
    }

    public RuntimeEOFException(final String s, final Throwable throwable) {
        super(s, throwable);
    }

    public RuntimeEOFException(final Throwable throwable) {
        super(throwable);
    }
}
