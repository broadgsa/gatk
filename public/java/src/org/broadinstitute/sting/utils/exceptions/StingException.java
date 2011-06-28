package org.broadinstitute.sting.utils.exceptions;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 8:11:12 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 6, 2009
 * <p/>
 * Class StingException
 * <p/>
 * This exception allows us to filter out exceptions that come from core GATK code, and those that
 * are not homegrown..
 */
public class StingException extends RuntimeException {
    public StingException(String msg) {
        super(msg);
    }

    public StingException(String message, Throwable throwable) {
        super(message, throwable);
    }
}

