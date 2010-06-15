package org.broadinstitute.sting.queue.util

import org.apache.log4j._

/**
 * A mixin to add logging to a class
 */
trait Logging {
  private val className = this.getClass.getName
  protected lazy val logger = configuredLogger

  def configuredLogger = {
    Logging.configureLogging
    Logger.getLogger(className)
  }
}

object Logging {
  private var configured = false
  private var isDebug = false
  def configureLogging = {
    if (!configured) {
      var root = Logger.getRootLogger
      root.addAppender(new ConsoleAppender(new PatternLayout("%-5p %d{HH:mm:ss,SSS} - %m %n")))
      root.setLevel(if(isDebug) Level.DEBUG else Level.INFO)
      configured = true
    }
  }

  def enableDebug = {isDebug = true; Logger.getRootLogger.setLevel(Level.DEBUG)}
}
