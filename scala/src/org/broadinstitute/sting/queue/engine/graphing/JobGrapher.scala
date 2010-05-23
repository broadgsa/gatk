package org.broadinstitute.sting.queue.engine.graphing

import org.broadinstitute.sting.queue.QException
import org.jgrapht.graph.SimpleDirectedGraph
import org.jgrapht.alg.BellmanFordShortestPath
import collection.JavaConversions._
import org.broadinstitute.sting.queue.engine.{QCommand, QFile, QRule}
import org.broadinstitute.sting.queue.engine.scheduling.{ExecEdge, ResourceEdge, ResourceNode}
import org.broadinstitute.sting.queue.util.Logging

/**
 * A basic job grapher.
 * Limitiations:
 *  - Only walks along graphs with rules that have a single input and a single output.
 */
class JobGrapher(
    private val inputFiles: List[String],
    private val argMap: Map[String, String],
    private val rules: List[QRule],
    private val sourceFiles: List[QFile],
    private val targetFiles: List[QFile]) extends Logging {

  private val modelGraph = new SimpleDirectedGraph[QFile, QCommand](classOf[QCommand])
  private val jobGraph = new SimpleDirectedGraph[ResourceNode, ResourceEdge](classOf[ResourceEdge])

  createModelGraph()
  createJobGraph()

  def name = this.getClass.getName
  def jobs = jobGraph

  private def createJobGraph() = {
    var missingPaths = List.empty[Tuple2[QFile,QFile]]
    for (sourceFile <- sourceFiles) {
      for (targetFile <- targetFiles) {
        var shortestPath = BellmanFordShortestPath.findPathBetween(modelGraph, sourceFile, targetFile)
        if (shortestPath == null)
          missingPaths = missingPaths ::: List((sourceFile, targetFile))
        else
          addPaths(shortestPath)
      }
    }

    for ((sourceFile, targetFile) <- missingPaths) {
      logger.error(String.format("No command path found between %s --> %s", sourceFile, targetFile))
    }

    if (missingPaths.size > 0)
      throw new QException("Not all inputs and outputs found in the pipeline graph")
  }

  private def createModelGraph() = {
    for (rule <- rules) {
      if (rule.inputs.size != 1 || (rule.outputs.size != 1))
        throw new QException(this.name + " can only process rules with a single input and a single output.  " +
                "inputs: " + rule.inputs + ", outputs: " + rule.outputs + ", command: " + rule.command)
      var source = rule.inputs.head
      var target = rule.outputs.head
      modelGraph.addVertex(source)
      modelGraph.addVertex(target)
      modelGraph.addEdge(source, target, rule.command)
    }
  }

  private def addPaths(shortestPath: java.util.List[QCommand]) {
    for (inputFile <- inputFiles)
      if (modelGraph.getEdgeSource(shortestPath.head).matchesFile(inputFile))
        addPath(shortestPath, inputFile)
  }

  private def addPath(shortestPath: java.util.List[QCommand], inputFile: String) = {
    var sourceFile = inputFile
    for (command <- shortestPath) {
      val source = modelGraph.getEdgeSource(command)
      val target = modelGraph.getEdgeTarget(command)
      val baseName = source.baseName(sourceFile)
      val targetFile = target.fullName(baseName)
      val resourceSource = new ResourceNode(Map(source.extension -> sourceFile))
      val resourceTarget = new ResourceNode(Map(target.extension -> targetFile))
      val resourceEdge = new ExecEdge(argMap, command)
      jobGraph.addVertex(resourceSource)
      jobGraph.addVertex(resourceTarget)
      jobGraph.addEdge(resourceSource, resourceTarget, resourceEdge)
      sourceFile = targetFile
    }
  }
}
