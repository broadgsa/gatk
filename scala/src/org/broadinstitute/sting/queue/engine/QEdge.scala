package org.broadinstitute.sting.queue.engine

import java.io.File

/**
 * An edge in the QGraph
 */
trait QEdge {
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
