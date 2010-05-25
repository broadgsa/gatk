package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.DirectedGraph
import org.broadinstitute.sting.queue.util.Logging
import collection.JavaConversions._
import org.broadinstitute.sting.queue.QArguments

abstract class JobScheduler(protected val jobGraph: DirectedGraph[ResourceNode, ResourceEdge],
        protected val qArgs: QArguments) extends Logging {

  def runJobs
  def numJobs = jobGraph.edgeSet.size

  /**
   * Emulates storing of properties per node by looking up values on
   * the current edge/target-node or any preceding nodes in the graph.
   */
  def lookup(edge: ResourceEdge, key: String, default: String) : String = {
    lookupRecursive(edge, key) match {
      case Some(value) => value
      case None => qArgs.argMap.getOrElse(key, default)
    }
  }

  private def lookupRecursive(edge: ResourceEdge, key: String) : Option[String] = {
    var value = edge.lookup(key)
    if (value.isDefined)
      return value

    value = this.jobGraph.getEdgeTarget(edge).lookup(key)
    if (value.isDefined)
      return value

    return lookupRecursive(this.jobGraph.getEdgeSource(edge), key)
  }

  private def lookupRecursive(node: ResourceNode, key: String) : Option[String] = {
    var value = node.lookup(key)
    if (value.isDefined)
      return value

    for (edge <- this.jobGraph.incomingEdgesOf(node)) {
      value = lookupRecursive(edge, key)
      if (value.isDefined)
        return value
    }

    return None
  }
}
