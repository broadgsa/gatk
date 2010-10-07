package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{ArgumentSource, Input, Output}
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base class for Gather command line functions.
 */
trait GatherFunction extends QFunction {
  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[File] = Nil

  @Output(doc="The original output of the scattered function")
  var originalOutput: File = _

  /**
   * Sets the original function used to create this scatter function.
   * @param originalFunction The ScatterGatherableFunction.
   * @param gatherField The field being gathered.
   */
  def setOriginalFunction(originalFunction: ScatterGatherableFunction, gatherField: ArgumentSource) = {}

  /**
   * Sets the clone function creating one of the inputs for this gather function.
   * @param cloneFunction The clone wrapper for the original ScatterGatherableFunction.
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   * @param gatherField The field to be gathered.
   */
  def setCloneFunction(cloneFunction: CloneFunction, index: Int, gatherField: ArgumentSource) = {}
}
