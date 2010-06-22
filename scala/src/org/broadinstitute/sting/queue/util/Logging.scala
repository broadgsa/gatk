package org.broadinstitute.sting.queue.util

import org.apache.log4j._

/**
 * A mixin to add logging to a class
 */
trait Logging {
  private val className = this.getClass.getName
  protected lazy val logger = {
    Logging.configureLogging
    Logger.getLogger(className)
  }
}

object Logging {
  private var configured = false
  private var level = Level.INFO
  def configureLogging = {
    if (!configured) {
      var root = Logger.getRootLogger
      root.addAppender(new ConsoleAppender(new PatternLayout("%-5p %d{HH:mm:ss,SSS} - %m %n")))
      root.setLevel(level)
      configured = true
    }
  }

  def setDebug = setLevel(Level.DEBUG)
  def setTrace = setLevel(Level.TRACE)
  private def setLevel(level: Level) = {this.level = level; Logger.getRootLogger.setLevel(level)}
}
