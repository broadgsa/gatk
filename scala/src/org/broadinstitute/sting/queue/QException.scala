package org.broadinstitute.sting.queue

import org.broadinstitute.sting.utils.exceptions.StingException

class QException(private val message: String, private val throwable: Throwable = null)
  extends StingException(message, throwable)
