package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.commandline.{Input, Output}
import java.io.File

trait ScatterFunction extends CommandLineFunction {
  type ScatterType

  @Input(doc="Original input to scatter")
  var originalInput: ScatterType = _

  @Input(doc="Temporary directories for each scatter part")
  var tempDirectories: List[File] = Nil

  @Output(doc="Scattered parts of the original input, one per temp directory")
  var scatterParts: List[ScatterType] = Nil

  def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {}
}
