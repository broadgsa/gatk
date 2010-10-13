package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.QException

/**
 * Thrown after giving up on retrying.
 */
class RetryException(private val message: String, private val throwable: Throwable)
        extends QException(message, throwable) {}
