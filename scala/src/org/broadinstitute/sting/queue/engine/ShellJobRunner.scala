package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.{JobExitException, Logging, ShellJob}
import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Runs jobs one at a time locally
 */
class ShellJobRunner(function: CommandLineFunction) extends JobRunner with Logging {
  private var runStatus: RunnerStatus.Value = _

  /**
   * Runs the function on the local shell.
   * @param function Command to run.
   */
  def start() = {
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
    if (function.jobErrorFile != null) {
      logger.info("Errors written to " + function.jobErrorFile)
    } else {
      if (logger.isDebugEnabled)
        logger.info("Errors also written to " + function.jobOutputFile)
    }

    function.doneOutputs.foreach(_.delete)
    function.failOutputs.foreach(_.delete)
    runStatus = RunnerStatus.RUNNING
    try {
        job.run()
        function.doneOutputs.foreach(_.createNewFile)
        runStatus = RunnerStatus.DONE
        logger.info("Done: " + function.commandLine)
    } catch {
      case e: JobExitException =>
        runStatus = RunnerStatus.FAILED
        try {
          function.failOutputs.foreach(_.createNewFile)
        } catch {
          case _ => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + function.commandLine, e)
    }
  }

  def status = runStatus
}
