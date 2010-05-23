package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.DirectedGraph
import org.apache.commons.lang.text.{StrLookup, StrSubstitutor}
import collection.JavaConversions._
import org.broadinstitute.sting.queue.engine.QCommand

class ExecEdge(val args: Map[String, String], private val command: QCommand) extends ResourceEdge {
  private var convertedCommandString: String = _
  def commandString = convertedCommandString

  override def traverse(graph: DirectedGraph[ResourceNode, ResourceEdge]) = {
    // Lookup any variable using the target node, or any of it's input nodes.
    val sub = new StrSubstitutor(new NodeLookup(graph.getEdgeTarget(this), graph))
    convertedCommandString = sub.replace(command.commandString)
  }

  class NodeLookup(private val targetNode: ResourceNode, private val graph: DirectedGraph[ResourceNode, ResourceEdge]) extends StrLookup {

    def lookup(key: String) = {
      var value: String = null
      if (args.contains(key))
        value = args(key)
      else
        value = lookup(key, targetNode)
      value
    }

    private def lookup(key: String, node: ResourceNode): String = {
      var value: String = null
      if (node.resources.contains(key)) {
        value = node.resources(key)
      } else {
        for (edge <- graph.incomingEdgesOf(node)) {
          lookup(key, graph.getEdgeSource(edge)) match {
            case null => {}
            case found => value = found
          }
        }
      }
      value
    }
  }
}
