package org.broadinstitute.sting.queue.function

/**
 * Utility class to map a set of inputs to set of outputs.
 * The QGraph uses this function internally to map between user defined functions.
 */
class MappingFunction(val inputs: Set[Any], val outputs: Set[Any]) extends QFunction {
  override def toString = "<map>" // For debugging
}
