package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.{Logging, ShellJob}
import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Runs jobs one at a time locally
 */
class ShellJobRunner extends JobRunner with Logging {
  /**
   * Runs the function on the local shell.
   * @param function Command to run.
   * @param qGraph graph that holds the job, and if this is a dry run.
   */
  def run(function: CommandLineFunction, qGraph: QGraph) = {
    val job = new ShellJob
    job.command = function.commandLine
    job.workingDir = function.commandDirectory
    job.outputFile = function.jobOutputFile
    job.errorFile = function.jobErrorFile

    if (logger.isDebugEnabled) {
      logger.debug(function.commandDirectory + " > " + function.commandLine)
    } else {
      logger.info(function.commandLine)
    }
    logger.info("Output written to " + function.jobOutputFile)
    if (function.jobErrorFile != null) {
      logger.info("Errors written to " + function.jobErrorFile)
    } else {
      if (logger.isDebugEnabled)
        logger.info("Errors also written to " + function.jobOutputFile)
    }

    if (!qGraph.dryRun)
      job.run
  }
}
