package org.broadinstitute.sting.queue.engine

/**
 * Utility class to map a set of inputs to set of outputs.
 * The QGraph uses this function internally to map between user defined functions.
 */
class MappingEdge(val inputs: QNode, val outputs: QNode) extends QEdge {
  /**
   * For debugging purposes returns <map>.
   * @return <map>
   */
  override def toString = "<map>"
  override def dotString = "<map>"
}
