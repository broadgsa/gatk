package org.broadinstitute.sting.utils;
/**
 * User: hanna
 * Date: Jun 2, 2009
 * Time: 11:43:48 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Indicates that something was wrong with in the parameters passed to create a GenomeLoc...
 * bad sequence id out of bounds, etc.
 */

public class MalformedGenomeLocException extends StingException {
    /**
     * Create a new MalformedGenomeLocException with the given message.  Does not preserve the existing stack trace.
     * @param message The message.
     */
    public MalformedGenomeLocException( String message ) {
        super(message);
    }

    /**
     * Create a new MalformedGenomeLocException with the given message and root cause.
     * @param message The message.
     * @param t The root cause.
     */
    public MalformedGenomeLocException( String message, Throwable t ) {
        super(message,t);
    }
}
