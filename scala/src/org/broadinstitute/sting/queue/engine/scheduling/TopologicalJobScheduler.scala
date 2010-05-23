package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.graph.SimpleDirectedGraph
import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.event.{EdgeTraversalEvent, TraversalListenerAdapter}
import collection.JavaConversions._

/**
 * Loops over the job graph running jobs as the edges are traversed
 */
abstract class TopologicalJobScheduler(jobGraph: SimpleDirectedGraph[ResourceNode, ResourceEdge]) extends JobScheduler(jobGraph) {

  protected val iterator = new TopologicalOrderIterator(this.jobGraph)

  iterator.addTraversalListener(new TraversalListenerAdapter[ResourceNode, ResourceEdge] {
    override def edgeTraversed(event: EdgeTraversalEvent[ResourceNode, ResourceEdge]) = {
      traversed(event.getEdge)
    }
  })

  override def runJobs = {
    logger.info(String.format("Running %s jobs.", this.numJobs.toString))
    for (target <- iterator) {
      // Do nothing for now, let event handler respond
    }
  }

  protected def traversed(edge: ResourceEdge) = {
    edge.traverse(this.jobGraph)
    edge match {
      case exec: ExecEdge => traversedExec(exec)
    }
  }

  protected def traversedExec(exec: ExecEdge)
}
