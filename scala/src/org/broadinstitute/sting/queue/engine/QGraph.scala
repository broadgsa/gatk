package org.broadinstitute.sting.queue.engine

import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.immutable.ListMap
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.{MappingFunction, CommandLineFunction, QFunction}

class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  val jobGraph = new SimpleDirectedGraph[QNode, QFunction](classOf[QFunction])
  def numJobs = JavaConversions.asSet(jobGraph.edgeSet).filter(_.isInstanceOf[CommandLineFunction]).size

  def add(command: CommandLineFunction) {
    add(command, true)
  }

  /**
   * Looks through functions with multiple inputs and outputs and adds mapping functions for single inputs and outputs.
   */
  def fillIn = {
    // clone since edgeSet is backed by the graph
    for (function <- JavaConversions.asSet(jobGraph.edgeSet).clone) {
      val inputs = function.inputs
      val outputs = function.outputs

      if (inputs.size > 1)
        for ((name, input) <- inputs)
          addNullEdge(ListMap(name -> input), inputs)

      if (outputs.size > 1)
        for ((name, output) <- outputs)
          addNullEdge(outputs, ListMap(name -> output))
    }
  }

  def run = {
    var isReady = true
    for (function <- JavaConversions.asSet(jobGraph.edgeSet)) {
      val missingValues = function.missingValues
      if (missingValues.size > 0) {
        isReady = false
        logger.error(function match {
          case cmd: CommandLineFunction => "Missing values for function: %s".format(cmd.commandLine)
          case x => "Missing values:"
        })
        for (missing <- missingValues) {
          logger.error("  " + missing)
        }
      }
    }
    
    if (isReady || this.dryRun)
      (new TopologicalJobScheduler(this) with LsfJobRunner).runJobs
  }

  private def add(f: QFunction, replace: Boolean) {
    val inputs = QNode(f.inputs.values.filter(_ != null).toSet)
    val outputs = QNode(f.outputs.values.filter(_ != null).toSet)
    jobGraph.addVertex(inputs)
    jobGraph.addVertex(outputs)
    if (replace)
      jobGraph.removeAllEdges(inputs, outputs)
    jobGraph.addEdge(inputs, outputs, f)
  }

  private def addNullEdge(input: ListMap[String, Any], output: ListMap[String, Any]) = {
    add(new MappingFunction(input, output), false)
  }
}
