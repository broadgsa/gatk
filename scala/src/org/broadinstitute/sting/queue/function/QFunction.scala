package org.broadinstitute.sting.queue.function

import scala.collection.immutable.ListMap

/**
 * The base interface for all functions in Queue.
 * Inputs and outputs are specified as ListMaps of name -> value.
 * The names are used for debugging.
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
   * ListMap of name -> value inputs for this function.
   */
  def inputs: ListMap[String, Any]

  /**
   * ListMap of name -> value outputs for this function.
   */
  def outputs: ListMap[String, Any]
}
