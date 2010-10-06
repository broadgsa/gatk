package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.{QFunction, CommandLineFunction}

/**
 * Only logs the command to run.  Doesn't actually run it.
 */
class DryRunner(function: QFunction) extends JobRunner with Logging {
  /**
   * Dry runs the function logging the command lines.
   * @param function Command to run.
   */
  // TODO: Why do we need the dry runner?  Can we just use a dryRun() method to log per JobRunner?
  def start() = {
    function match {
      case clf: CommandLineFunction => {
        if (logger.isDebugEnabled) {
          logger.debug(clf.commandDirectory + " > " + clf.commandLine)
        } else {
          logger.info(clf.commandLine)
        }
        logger.info("Output written to " + clf.jobOutputFile)
        if (clf.jobErrorFile != null) {
          logger.info("Errors written to " + clf.jobErrorFile)
        } else {
          if (logger.isDebugEnabled)
            logger.info("Errors also written to " + clf.jobOutputFile)
        }
      }
      case qFunction => {
        logger.info(qFunction.description)
      }
    }
  }

  def status = RunnerStatus.DONE
}
