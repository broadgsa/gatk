package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.InProcessFunction
import java.util.Date
import org.broadinstitute.sting.utils.Utils
import org.apache.commons.io.{IOUtils, FileUtils}
import java.io.PrintStream

/**
 * Runs a function that executes in process and does not fork out an external process.
 */
class InProcessRunner(val function: InProcessFunction) extends JobRunner[InProcessFunction] {
  private var runStatus: RunnerStatus.Value = _

  def start() {
    getRunInfo.startTime = new Date()
    getRunInfo.exechosts = Utils.resolveHostname()
    runStatus = RunnerStatus.RUNNING

    function.jobOutputStream = new PrintStream(FileUtils.openOutputStream(function.jobOutputFile))
    function.jobErrorStream = {
      if (function.jobErrorFile != null)
        new PrintStream(FileUtils.openOutputStream(function.jobErrorFile))
      else
        function.jobOutputStream
    }
    try {
      function.run()
      function.jobOutputStream.println("%s%nDone.".format(function.description))
    } finally {
      IOUtils.closeQuietly(function.jobOutputStream)
      if (function.jobErrorFile != null)
        IOUtils.closeQuietly(function.jobErrorStream)
    }

    runStatus = RunnerStatus.DONE
    getRunInfo.doneTime = new Date()
  }

  def status = runStatus
}
