package org.broadinstitute.sting.queue.function

import java.io.File

/**
 * Utility class to map a set of inputs to set of outputs.
 * The QGraph uses this function internally to map between user defined functions.
 */
class MappingFunction(val inputs: Set[File], val outputs: Set[File]) extends QFunction {
  /**
   * For debugging purposes returns <map>.
   * @returns <map>
   */
  override def toString = "<map>"
}
