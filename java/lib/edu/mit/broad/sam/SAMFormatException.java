/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

/**
 * Thrown when a SAM file being read (text or binary) looks bad.
 */
public class SAMFormatException extends RuntimeException {
    public SAMFormatException() {
    }

    public SAMFormatException(final String s) {
        super(s);
    }

    public SAMFormatException(final String s, final Throwable throwable) {
        super(s, throwable);
    }

    public SAMFormatException(final Throwable throwable) {
        super(throwable);
    }
}
