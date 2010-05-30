package org.broadinstitute.sting.queue.engine.graphing

import org.jgrapht.graph.SimpleDirectedGraph
import org.jgrapht.alg.BellmanFordShortestPath
import org.broadinstitute.sting.queue.{QArguments, QException}
import collection.mutable.ListBuffer
import collection.JavaConversions._
import org.broadinstitute.sting.queue.engine.scheduling.{ResourceEdge, ExecEdge, MapResourceNode}
import org.broadinstitute.sting.queue.engine.{QModelEdge, QCommand, QFile, QRule}

/**
 * Converts a set of rules provided by the user and a list of files into a graph of jobs to run
 */
class TreeJobGrapher extends JobGrapher {
  private val modelGraph = new SimpleDirectedGraph[List[QFile], QModelEdge](classOf[QModelEdge])

  private val rules = new ListBuffer[QRule]
  private var sourceFiles = List.empty[QFile]
  private var targetFiles = List.empty[QFile]

  // Used to tag a model edge for Element <-> List
  private class QCollectionEdge extends QModelEdge

  override protected def createJobGraph() = {
    createModelGraph()
    
    var missingPaths = List.empty[Tuple2[QFile,QFile]]
    for (sourceFile <- sourceFiles) {
      for (targetFile <- targetFiles) {
        var shortestPath = BellmanFordShortestPath.findPathBetween(modelGraph, List(sourceFile), List(targetFile))
        if (shortestPath == null)
          missingPaths = missingPaths ::: List((sourceFile, targetFile))
        else
          addPaths(shortestPath, qArgs)
      }
    }

    for ((sourceFile, targetFile) <- missingPaths) {
      logger.error(String.format("No command path found between %s --> %s", sourceFile, targetFile))
    }

    if (missingPaths.size > 0)
      throw new QException("Not all inputs and outputs found in the pipeline graph")
  }

  private def createModelGraph() = {

    // Look for rules with more than one input or output and add
    // internal dependencies between the elements and the list
    for (rule <- rules) {
      if (rule.inputs.size > 1) {
        for (input <- rule.inputs) {
          modelGraph.addVertex(List(input))
          modelGraph.addVertex(rule.inputs)
          modelGraph.addEdge(List(input), rule.inputs, new QCollectionEdge)
        }
      }
      if (rule.outputs.size > 1) {
        for (output <- rule.outputs) {
          modelGraph.addVertex(rule.outputs)
          modelGraph.addVertex(List(output))
          modelGraph.addEdge(rule.outputs, List(output), new QCollectionEdge)
        }
      }
    }

    // Add the explicit rules
    for (rule <- rules) {
      modelGraph.addVertex(rule.inputs)
      modelGraph.addVertex(rule.outputs)
      modelGraph.addEdge(rule.inputs, rule.outputs, rule.command)
    }
  }

  private def addPaths(shortestPath: java.util.List[QModelEdge], qArgs: QArguments) {
    for (inputFile <- qArgs.inputPaths.map(_.getAbsolutePath)) {
      val source = modelGraph.getEdgeSource(shortestPath.head).head
      if (source.matchesFile(inputFile)) {
        val baseName = source.baseName(inputFile)
        val target = modelGraph.getEdgeTarget(shortestPath.last).head
        addPathsToTarget(baseName, List(target))
      }
    }
  }

  private def addPathsToTarget(baseName: String, targets: List[QFile]) : Unit = {
    for (command <- modelGraph.incomingEdgesOf(targets)) {
      val sources = modelGraph.getEdgeSource(command)
      addJobGraphEdge(baseName, sources, targets, command)
      addPathsToTarget(baseName, sources)
    }
  }

  private def addJobGraphEdge(baseName: String, sources: List[QFile], targets: List[QFile], command: QModelEdge) {
    val resourceSource = new MapResourceNode(mapFiles(baseName, sources))
    val resourceTarget = new MapResourceNode(mapFiles(baseName, targets))
    val resourceEdge = command match {
      case qCommand: QCommand => new ExecEdge(qCommand.commandString)
      case qTransition: QCollectionEdge => new ResourceEdge
    }
    jobGraph.addVertex(resourceSource)
    jobGraph.addVertex(resourceTarget)
    jobGraph.addEdge(resourceSource, resourceTarget, resourceEdge)
  }

  /**
   * Creates a mapping of the files based on the baseName.
   *   key: file extension
   *   value: full name
   * Used by the JobScheduler to lookup values as the commands are expanded at exec time.
   */
  private def mapFiles(baseName: String, files: List[QFile]) : Map[String, String] = {
    Map(files.map(file => (file.extension, file.fullName(baseName))):_*)
  }
}

/**
 * Syntactic sugar for filling in a pipeline using a Scala script.
 */
object TreeJobGrapher {
  private val grapher = new TreeJobGrapher

  /**
   * Sugar that allows addRule( inputs -> outputs, command )
   */
  def addRule(rule: (Any, Any), commandString: String): Unit = {
    val inputs = QFile.getFiles(rule._1)
    val outputs = QFile.getFiles(rule._2)
    val command = new QCommand(commandString)
    addRule(inputs, outputs, command)
  }

  private def addRule(inputs: List[QFile], outputs: List[QFile], command: QCommand): Unit = {
    grapher.rules += new QRule(inputs, outputs, command)
  }

  def run(args: Array[String], sources: Any, targets: Any) = {
    grapher.qArgs = new QArguments(args)
    grapher.sourceFiles = QFile.getFiles(sources)
    grapher.targetFiles = QFile.getFiles(targets)
    grapher.run()
  }
}
