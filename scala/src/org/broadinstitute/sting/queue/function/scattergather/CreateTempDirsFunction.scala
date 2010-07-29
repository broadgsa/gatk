package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.commandline.{Output, Input}

class CreateTempDirsFunction extends CommandLineFunction {
  @Input(doc="Original inputs to the scattered function")
  var originalInputs: Set[Any] = Set.empty[Any]

  @Output(doc="Temporary directories to create")
  var tempDirectories: List[File] = Nil

  def commandLine = "mkdir -pv%s".format(repeat(" '", tempDirectories, "'"))
}
