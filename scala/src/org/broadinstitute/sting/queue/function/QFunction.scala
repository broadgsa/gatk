package org.broadinstitute.sting.queue.function

import java.io.File

/**
 * The base interface for all functions in Queue.
 * Inputs and outputs are specified as Sets of values.
 * Inputs are matched to other outputs by using .equals()
 */
trait QFunction {
  /**
   * After a function is frozen no more updates are allowed by the user.
   * The function is allow to make necessary updates internally to make sure
   * the inputs and outputs will be equal to other inputs and outputs.
   */
  def freeze = {}

  /**
   * Set of inputs for this function.
   */
  def inputs: Set[File]

  /**
   * Set of outputs for this function.
   */
  def outputs: Set[File]

  /**
   * The function description in .dot files
   */
  def dotString = ""
}
