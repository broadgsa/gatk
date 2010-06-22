package org.broadinstitute.sting.queue.engine

import scala.collection.immutable.ListMap

/**
 * Represents a state between QFunctions the directed acyclic QGraph
 */
case class QNode (private val items: Set[Any]) {
  /**
   * Used during QGraph error reporting.
   * The EdgeFactory uses the valueMap to create new edges for the CycleDetector.
   */
  def valueMap = {
    var map = ListMap.empty[String, Any]
    for (item <- items)
      if (item != null)
        map += item.toString -> item
    map
  }
}
