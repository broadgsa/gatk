package org.broadinstitute.sting.queue.function

import java.io.File

/** An internal class that is used by bsub to wait on all other jobs before exiting. */
class DispatchWaitFunction extends CommandLineFunction {
  /**
   * Returns the command line "echo".
   * @return echo
   */
  def commandLine = "echo"

  jobQueue = "hour"
  jobOutputFile = File.createTempFile("Q-wait", ".out")
}
