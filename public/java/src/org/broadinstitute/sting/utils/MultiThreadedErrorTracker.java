package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * A utility to track exceptions that occur across threads.
 *
 * Uses a notify mechanism so that multiple threads can tell the tracker that an
 * error has occurred, and a master thread can monitor this object for an error
 * occurring and take appropriate action.  Only maintains the first
 * error to reach the tracker.
 *
 * Refactored from HierarchicalMicroScheduler
 *
 * User: depristo
 * Date: 9/19/12
 * Time: 11:20 AM
 */
public class MultiThreadedErrorTracker {
    /**
     * An exception that's occurred.  If null, no exception has occurred.
     */
    private RuntimeException error = null;

    /**
     * Convenience function to check, and throw, an error is one is pending
     */
    public synchronized void throwErrorIfPending() {
        if (hasAnErrorOccurred())
            throw getError();
    }

    /**
     * Detects whether an execution error has occurred.
     * @return True if an error has occurred.  False otherwise.
     */
    public synchronized boolean hasAnErrorOccurred() {
        return error != null;
    }

    /**
     * Retrieve the error that has occurred.
     *
     * @throws ReviewedStingException if no error has occurred.
     * @return
     */
    public synchronized RuntimeException getError() {
        if(!hasAnErrorOccurred())
            throw new ReviewedStingException("User has attempted to retrieve a traversal error when none exists");
        return error;
    }

    /**
     * Notify this error tracker that an error has occurs.  Only updates the tracked
     * error if it is currently null (i.e., no error has been already reported).  So
     * calling this successively with multiple errors only keeps the first, which is the
     * right thing to do as the initial failure is usually the meaningful one, but
     * generates a cascade of failures as other subsystems fail.
     */
    public synchronized RuntimeException notifyOfError(Throwable error) {
        if ( this.error == null )
            this.error = toRuntimeException(error);

        return this.error;
    }

    /**
     * Convert error to a Runtime exception, or keep as is if it already is one
     *
     * @param error the error that has occurred
     * @return the potentially converted error
     */
    private RuntimeException toRuntimeException(final Throwable error) {
        // If the error is already a Runtime, pass it along as is.  Otherwise, wrap it.
        if (error instanceof RuntimeException)
            return (RuntimeException)error;
        else
            return new ReviewedStingException("An error occurred during the traversal.  Message=" + error.getMessage(), error);
    }
}
