package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.IOUtils
import java.io.{PrintWriter, StringWriter}
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base interface for job runners.
 */
trait JobRunner[TFunction <: QFunction] {
  /**
   * Runs the function.
   * After the function returns the status of the function should
   * be RUNNING, FAILED, or DONE.
   * @param function Command to run.
   */
  def start()

  /**
   * Returns the current run status.
   * Must only be called AFTER start().
   * @return RUNNING, DONE, or FAILED.
   */
  def status: RunnerStatus.Value

  /**
   * Returns the function to be run.
   */
  def function: TFunction

  protected def writeDone() = {
    val content = "%s%nDone.".format(function.description)
    IOUtils.writeContents(function.jobOutputFile, content)
  }

  protected def writeError(content: String) = {
    IOUtils.writeContents(functionErrorFile, content)
  }

  /**
   * Writes the stack trace to the error file.
   */
  protected def writeStackTrace(e: Throwable) = {
    val stackTrace = new StringWriter
    val printWriter = new PrintWriter(stackTrace)
    printWriter.println(function.description)
    e.printStackTrace(printWriter)
    printWriter.close
    IOUtils.writeContents(functionErrorFile, stackTrace.toString)
  }

  /**
   * Calls back to a hook that an expert user can setup to modify a job.
   * @param value Value to modify.
   */
  protected def updateJobRun(value: Any) = {
    val updater = function.updateJobRun
    if (updater != null)
      if (updater.isDefinedAt(value))
        updater(value)
  }

  private def functionErrorFile = if (function.jobErrorFile != null) function.jobErrorFile else function.jobOutputFile
}
