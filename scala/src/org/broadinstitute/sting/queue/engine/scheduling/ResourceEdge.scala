package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.DirectedGraph

abstract class ResourceEdge {
  def traverse(graph: DirectedGraph[ResourceNode, ResourceEdge]): Unit
}
