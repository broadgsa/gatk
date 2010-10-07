package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{ArgumentSource, Input, Output}
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base class for Scatter command line functions.
 */
trait ScatterFunction extends QFunction {
  @Input(doc="Original input to scatter")
  var originalInput: File = _

  @Output(doc="Scattered parts of the original input, one per temp directory")
  var scatterParts: List[File] = Nil

  @Input(doc="Temporary directories for each scatter part")
  var tempDirectories: List[File] = Nil

  /**
   * Sets the original function used to create this scatter function.
   * @param originalFunction The ScatterGatherableFunction.
   * @param scatterField The field being scattered.
   */
  def setOriginalFunction(originalFunction: ScatterGatherableFunction, scatterField: ArgumentSource) = {}

  /**
   * Sets the clone function using one of the outputs of this scatter function.
   * @param cloneFunction The clone wrapper for the original ScatterGatherableFunction.
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   * @param scatterField The field being scattered.
   */
  def setCloneFunction(cloneFunction: CloneFunction, index: Int, scatterField: ArgumentSource) = {}
}
