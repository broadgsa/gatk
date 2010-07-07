package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.{CommandLineFunction}
import org.broadinstitute.sting.commandline.{Input, Output}

/**
 * Base class for Gather command line functions.
 * NOTE: Using an abstract class instead of a trait due to scala parameterized type erasure on traits.
 */
abstract class GatherFunction extends CommandLineFunction {
  type GatherType

  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[GatherType] = Nil

  @Output(doc="The original output of the scattered function")
  var originalOutput: GatherType = _

  def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {}
}
