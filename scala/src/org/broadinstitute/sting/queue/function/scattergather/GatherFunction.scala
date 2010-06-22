package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.{CommandLineFunction}
import org.broadinstitute.sting.queue.util.{Input, Output}

trait GatherFunction extends CommandLineFunction {
  type GatherType

  @Input
  var gatherParts: List[GatherType] = Nil

  @Output
  var originalOutput: GatherType = _

  def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {}
}
