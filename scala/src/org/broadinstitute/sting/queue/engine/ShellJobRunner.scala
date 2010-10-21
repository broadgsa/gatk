package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.{Logging, ShellJob}
import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Runs jobs one at a time locally
 */
class ShellJobRunner(val function: CommandLineFunction) extends JobRunner with Logging {
  private var runStatus: RunnerStatus.Value = _

  /**
   * Runs the function on the local shell.
   * @param function Command to run.
   */
  def start() = {
    try {
      val job = new ShellJob
      job.command = function.commandLine
      job.workingDir = function.commandDirectory
      job.outputFile = function.jobOutputFile
      job.errorFile = function.jobErrorFile

      if (logger.isDebugEnabled) {
        logger.debug("Starting: " + function.commandDirectory + " > " + function.commandLine)
      } else {
        logger.info("Starting: " + function.commandLine)
      }

      logger.info("Output written to " + function.jobOutputFile)
      if (function.jobErrorFile != null)
        logger.info("Errors written to " + function.jobErrorFile)

      function.deleteLogs()
      function.deleteOutputs()
      runStatus = RunnerStatus.RUNNING
      function.mkOutputDirectories()
      job.run()
      function.doneOutputs.foreach(_.createNewFile())
      runStatus = RunnerStatus.DONE
      logger.info("Done: " + function.commandLine)
    } catch {
      case e =>
        runStatus = RunnerStatus.FAILED
        try {
          function.failOutputs.foreach(_.createNewFile())
          writeStackTrace(e)
        } catch {
          case _ => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + function.commandLine, e)
    }
  }

  def status = runStatus
}
