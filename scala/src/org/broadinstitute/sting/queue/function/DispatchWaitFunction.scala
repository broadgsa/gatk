package org.broadinstitute.sting.queue.function

import java.io.File

class DispatchWaitFunction extends CommandLineFunction {
  def commandLine = "echo"

  jobQueue = "short"
  jobOutputFile = File.createTempFile("Q-wait", ".out")
  jobErrorFile = File.createTempFile("Q-wait", ".err")
}
