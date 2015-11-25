/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
     * @throws ReviewedGATKException if no error has occurred.
     * @return
     */
    public synchronized RuntimeException getError() {
        if(!hasAnErrorOccurred())
            throw new ReviewedGATKException("User has attempted to retrieve a traversal error when none exists");
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
            return new ReviewedGATKException("An error occurred during the traversal.  Message=" + error.getMessage(), error);
    }
}
