package org.broadinstitute.sting.queue.function

import java.io.File
import org.broadinstitute.sting.queue.util.Internal

trait DispatchFunction extends QFunction with MemoryLimitedFunction {
  def commandLine: String
  var commandDirectory: File

  var jobName: String = _
  var jobOutputFile: File = _
  var jobErrorFile: File = _

  @Internal
  var jobProject = "Queue"

  @Internal
  var jobQueue = "broad"
}
