package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.DirectedGraph
import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.event.{EdgeTraversalEvent, TraversalListenerAdapter}
import collection.JavaConversions._
import org.broadinstitute.sting.queue.QArguments

/**
 * Loops over the job graph running jobs as the edges are traversed
 */
abstract class TopologicalJobScheduler(jobGraph: DirectedGraph[ResourceNode, ResourceEdge], qArgs: QArguments)
        extends JobScheduler(jobGraph, qArgs) {

  protected val iterator = new TopologicalOrderIterator(this.jobGraph)

  iterator.addTraversalListener(new TraversalListenerAdapter[ResourceNode, ResourceEdge] {
    override def edgeTraversed(event: EdgeTraversalEvent[ResourceNode, ResourceEdge]) =
      event.getEdge.traverse(TopologicalJobScheduler.this)
  })

  override def runJobs = {
    logger.info(String.format("Running %s jobs.", this.numJobs.toString))
    for (target <- iterator) {
      // Do nothing for now, let event handler respond
    }
    logMissingKeys
  }
}
