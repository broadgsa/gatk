package org.broadinstitute.sting.queue.engine

/**
 * An edge in the QGraph
 */
trait QEdge {
  /**
   * List of inputs for this function sorted by path.
   */
  def inputs: QNode

  /**
   * List of outputs for this function sorted by path.
   */
  def outputs: QNode

  /**
   * The function description in .dot files
   */
  def dotString = ""

  override def hashCode = inputs.hashCode + outputs.hashCode

  override def equals(obj: Any) = {
    obj match {
      case other: QEdge =>
        this.inputs == other.inputs &&
        this.outputs == other.outputs
      case _ => false
    }
  }
}
