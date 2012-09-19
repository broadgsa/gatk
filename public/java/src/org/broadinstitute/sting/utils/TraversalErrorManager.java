package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 9/19/12
 * Time: 11:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraversalErrorManager {
    /**
     * An exception that's occurred in this traversal.  If null, no exception has occurred.
     */
    private RuntimeException error = null;

    public synchronized void throwErrorIfPending() {
        if (hasTraversalErrorOccurred())
            throw getTraversalError();
    }

    /**
     * Detects whether an execution error has occurred.
     * @return True if an error has occurred.  False otherwise.
     */
    public synchronized boolean hasTraversalErrorOccurred() {
        return error != null;
    }

    public synchronized RuntimeException getTraversalError() {
        if(!hasTraversalErrorOccurred())
            throw new ReviewedStingException("User has attempted to retrieve a traversal error when none exists");
        return error;
    }

    /**
     * Allows other threads to notify of an error during traversal.
     */
    public synchronized RuntimeException notifyOfTraversalError(Throwable error) {
        // If the error is already a Runtime, pass it along as is.  Otherwise, wrap it.
        this.error = toRuntimeException(error);
        return this.error;
    }

    private RuntimeException toRuntimeException(final Throwable error) {
        // If the error is already a Runtime, pass it along as is.  Otherwise, wrap it.
        if (error instanceof RuntimeException)
            return (RuntimeException)error;
        else
            return new ReviewedStingException("An error occurred during the traversal.  Message=" + error.getMessage(), error);
    }
}
