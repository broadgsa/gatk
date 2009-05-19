package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

/**
 * Generic class for handling misc parsing exceptions.
 */
public class ParseException extends StingException {
    public ParseException( String message ) {
        super( message );
    }
}

