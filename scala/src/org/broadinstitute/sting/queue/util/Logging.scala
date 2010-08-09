package org.broadinstitute.sting.queue.util

import org.apache.log4j._

/**
 * A mixin to add logging to a class
 */
trait Logging {
  private val className = this.getClass.getName
  protected lazy val logger = Logger.getLogger(className)
}
