package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.queue.util.IOUtils

/**
 * Runs a function that executes in process and does not fork out an external process.
 */
class InProcessRunner(val function: InProcessFunction) extends JobRunner[InProcessFunction] {
  private var runStatus: RunnerStatus.Value = _

  def start() = {
    runStatus = RunnerStatus.RUNNING
    function.run()
    val content = "%s%nDone.".format(function.description)
    IOUtils.writeContents(function.jobOutputFile, content)
    runStatus = RunnerStatus.DONE
  }

  def status = runStatus
}
