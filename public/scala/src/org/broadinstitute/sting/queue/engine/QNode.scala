package org.broadinstitute.sting.queue.engine

import java.io.File

/**
 * Represents a state between QFunctions the directed acyclic QGraph
 * @param files The list of files that represent this node state ordered by file name.
 */
class QNode (val id: Int, val files: List[File]) {
  override def equals(obj: Any) = {
    obj match {
      case other: QNode => this.id == other.id
      case _ => false
    }
  }

  override def hashCode = id

  override def toString = files.toString
}
