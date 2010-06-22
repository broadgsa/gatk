package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.Input
import java.io.File

class CleanupTempDirsFunction extends CommandLineFunction {
  @Input
  var originalOutputs: List[Any] = Nil

  @Input
  var tempDirectories: List[File] = Nil

  def commandLine = "rm -rf%s".format(repeat(" '", tempDirectories, "'"))
}
