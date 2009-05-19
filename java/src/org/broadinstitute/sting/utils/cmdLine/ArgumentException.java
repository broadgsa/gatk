package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

/**
 * Generic class for handling misc parsing exceptions.
 */
public class ArgumentException extends StingException {
    public ArgumentException( String message ) {
        super( message );
    }
}

