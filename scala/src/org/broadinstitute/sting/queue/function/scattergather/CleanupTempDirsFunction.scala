package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.commandline.Input
import java.io.File

class CleanupTempDirsFunction extends CommandLineFunction {
  @Input(doc="Original outputs of the gather functions")
  var originalOutputs: Set[Any] = Set.empty[Any]

  @Input(doc="Temporary directories to be deleted")
  var tempDirectories: List[File] = Nil

  def commandLine = "rm -rf%s".format(repeat(" '", tempDirectories, "'"))
}
