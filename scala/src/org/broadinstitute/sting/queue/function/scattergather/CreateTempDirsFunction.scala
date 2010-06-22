package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.{Output, Input}

class CreateTempDirsFunction extends CommandLineFunction {
  @Input
  var originalInputs: List[Any] = Nil

  @Output
  var tempDirectories: List[File] = Nil

  def commandLine = "mkdir%s".format(repeat(" '", tempDirectories, "'"))
}
