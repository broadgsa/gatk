package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.queue.util.Logging

/**
 * Runs a function that executes in process and does not fork out an external process.
 */
class InProcessRunner(function: InProcessFunction) extends JobRunner with Logging {
  private var runStatus: RunnerStatus.Value = _

  def start() = {
    if (logger.isDebugEnabled) {
      logger.debug("Starting: " + function.commandDirectory + " > " + function.description)
    } else {
      logger.info("Starting: " + function.description)
    }

    function.doneOutputs.foreach(_.delete)
    function.failOutputs.foreach(_.delete)
    runStatus = RunnerStatus.RUNNING
    try {
      function.run()
      function.doneOutputs.foreach(_.createNewFile)
      runStatus = RunnerStatus.DONE
      logger.info("Done: " + function.description)
    } catch {
      case e => {
        runStatus = RunnerStatus.FAILED
        try {
          function.failOutputs.foreach(_.createNewFile)
        } catch {
          case _ => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + function.description, e)
      }
    }
  }

  def status = runStatus
}
