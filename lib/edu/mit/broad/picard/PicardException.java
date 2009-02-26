/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard;

/**
 * Basic Picard runtime exception that, for now, does nothing much
 *
 * @author Kathleen Tibbetts
 */
public class PicardException extends RuntimeException
{
    public PicardException(String message) {
        super(message);
    }

    public PicardException(String message, Throwable throwable) {
        super(message, throwable);
    }

}
