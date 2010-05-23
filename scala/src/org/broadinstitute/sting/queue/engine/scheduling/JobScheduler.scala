package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.graph.SimpleDirectedGraph
import org.broadinstitute.sting.queue.util.Logging

abstract class JobScheduler(protected val jobGraph: SimpleDirectedGraph[ResourceNode, ResourceEdge])
        extends Logging {
  def runJobs
  def numJobs = jobGraph.edgeSet.size
}
