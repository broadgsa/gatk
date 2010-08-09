package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.event.{EdgeTraversalEvent, TraversalListenerAdapter}
import collection.JavaConversions._
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function._

/**
 * Loops over the job graph running jobs as the edges are traversed.
 * @param val The graph that contains the jobs to be run.
 */
abstract class TopologicalJobScheduler(private val qGraph: QGraph)
    extends ShellJobRunner with DispatchJobRunner with Logging {

  protected val iterator = new TopologicalOrderIterator(qGraph.jobGraph)

  iterator.addTraversalListener(new TraversalListenerAdapter[QNode, QFunction] {
    /**
     * As each edge is traversed, either dispatch the job or run it locally.
     * @param event Event holding the edge that was passed.
     */
    override def edgeTraversed(event: EdgeTraversalEvent[QNode, QFunction]) = event.getEdge match {
      case f: CommandLineFunction if (qGraph.bsubAllJobs) => dispatch(f, qGraph)
      case f: CommandLineFunction => run(f, qGraph)
      case f: MappingFunction => /* do nothing for mapping functions */
    }
  })

  /**
   * Runs the jobs by traversing the graph.
   */
  def runJobs = {
    logger.info("Number of jobs: %s".format(qGraph.numJobs))
    if (logger.isTraceEnabled)
      logger.trace("Number of nodes: %s".format(qGraph.jobGraph.vertexSet.size))
    var numNodes = 0
    for (target <- iterator) {
      if (logger.isTraceEnabled)
        logger.trace("Visiting: " + target)
      numNodes += 1
      // Do nothing for now, let event handler respond
    }
    if (logger.isTraceEnabled)
      logger.trace("Done walking %s nodes.".format(numNodes))

    if (qGraph.bsubAllJobs && qGraph.bsubWaitJobs) {
      logger.info("Waiting for jobs to complete.")
      val wait = new DispatchWaitFunction
      wait.freeze
      dispatch(wait, qGraph)
    }
  }
}
