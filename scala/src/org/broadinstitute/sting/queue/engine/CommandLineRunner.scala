package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.{Logging, ProcessUtils}
import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Runs jobs one at a time locally
 */
trait CommandLineRunner extends Logging {
  def run(function: CommandLineFunction, qGraph: QGraph) = {
    var commandLine = function.commandLine
    logger.info(commandLine)

    if (!qGraph.dryRun)
      ProcessUtils.runCommandAndWait(function.commandLine, function.commandDirectory)
  }
}
