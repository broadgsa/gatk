package org.broadinstitute.sting.queue.engine.graphing

import org.broadinstitute.sting.queue.engine.scheduling.{ResourceNode, ExecEdge}
import org.broadinstitute.sting.queue.engine.QCommand

class ExplicitJobGrapher extends JobGrapher

object ExplicitJobGrapher {
  private val grapher = new ExplicitJobGrapher

  def node() = new ResourceNode

  def addEdge(rule: (ResourceNode, ResourceNode), commandString: String): Unit = {
    addEdge(rule._1, rule._2, new QCommand(commandString))
  }

  private def addEdge(source: ResourceNode, target: ResourceNode, command: QCommand) = {
    val resourceEdge = new ExecEdge(command)
    grapher.jobGraph.addVertex(source)
    grapher.jobGraph.addVertex(target)
    grapher.jobGraph.addEdge(source, target, resourceEdge)
  }

  def run() = {
    grapher.run()
  }
}
