package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{Input, Output}
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base class for Scatter functions.
 */
trait ScatterFunction extends QFunction {
  @Input(doc="Original inputs to scatter")
  var originalInputs: Set[File] = _

  @Output(doc="Scattered parts of the original inputs, one set per temp directory")
  var scatterParts: List[File] = Nil

  /**
   * Returns true if the scatter function can scatter this original function.
   * @param originalFunction The original function to check.
   * @return true if the scatter function can scatter this original function.
   */
  def isScatterGatherable(originalFunction: ScatterGatherableFunction): Boolean

  /**
   * Sets the original ScatterGatherableFunction to be scattered.
   * @param originalFunction The original function to with inputs bind to this scatter function.
   */
  def setScatterGatherable(originalFunction: ScatterGatherableFunction)

  /**
   * Initializes the input fields for the clone function.
   * The input values should be set to their defaults
   * and may be changed by the user.
   * @param cloneFunction CloneFunction to initialize.
   * @param index The one based scatter index.
   */
  def initCloneInputs(cloneFunction: CloneFunction, index: Int)

  /**
   * Binds the input fields for the clone function to this scatter function.
   * The input values should be set to their absolute values and added
   * to scatter parts.
   * @param cloneFunction CloneFunction to bind.
   * @param index The one based scatter index.
   */
  def bindCloneInputs(cloneFunction: CloneFunction, index: Int)
}
