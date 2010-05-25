package org.broadinstitute.sting.queue.engine.graphing

import org.broadinstitute.sting.queue.QArguments
import org.broadinstitute.sting.queue.engine.scheduling._
import org.broadinstitute.sting.queue.util.Logging
import org.jgrapht.graph.SimpleDirectedGraph

abstract class JobGrapher() extends Logging {
  /**
   * Jobs to be run.
   * Can be populated adhoc or during createJobGraph()
   */
  protected val jobGraph = new SimpleDirectedGraph[ResourceNode, ResourceEdge](classOf[ResourceEdge])

  var qArgs: QArguments = _

  def run() = {
    createJobGraph()
    val scheduler = createScheduler()
    scheduler.runJobs
  }

  protected def createJobGraph() = {}

  private def createScheduler() : JobScheduler = {
    qArgs.useBsub match {
      case false => new SimpleJobScheduler(jobGraph, qArgs)
      case true => new DispatchJobScheduler(jobGraph, qArgs)
    }
  }
}
