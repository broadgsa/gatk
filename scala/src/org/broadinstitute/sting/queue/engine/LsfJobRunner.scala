package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util._

/**
 * Runs jobs on an LSF compute cluster.
 */
abstract class LsfJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {
  protected var runStatus: RunnerStatus.Value = _

  /** Job Id of the currently executing job. */
  var jobId = -1L

  // TODO: Full bsub command for debugging.
  protected def bsubCommand = "bsub " + function.commandLine
}
