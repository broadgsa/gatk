package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.{CommandLineFunction}
import org.broadinstitute.sting.commandline.{Input, Output}

trait GatherFunction extends CommandLineFunction {
  type GatherType

  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[GatherType] = Nil

  @Output(doc="The original output of the scattered function")
  var originalOutput: GatherType = _

  def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {}
}
