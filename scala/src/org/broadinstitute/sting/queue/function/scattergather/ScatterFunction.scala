package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.{Input, Output}
import java.io.File

trait ScatterFunction extends CommandLineFunction {
  type ScatterType

  @Input
  var originalInput: ScatterType = _

  @Input
  var tempDirectories: List[File] = Nil

  @Output
  var scatterParts: List[ScatterType] = Nil

  def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {}
}
