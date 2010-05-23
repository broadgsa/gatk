package org.broadinstitute.sting.queue

class QException(private val message: String, private val throwable: Throwable = null)
  extends RuntimeException(message, throwable)
