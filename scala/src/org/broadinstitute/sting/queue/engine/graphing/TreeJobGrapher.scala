package org.broadinstitute.sting.queue.engine.graphing

import org.jgrapht.graph.SimpleDirectedGraph
import org.broadinstitute.sting.queue.engine.{QCommand, QFile, QRule}
import org.broadinstitute.sting.queue.engine.scheduling.{ExecEdge, MapResourceNode, ResourceEdge, ResourceNode}
import org.jgrapht.alg.BellmanFordShortestPath
import org.broadinstitute.sting.queue.{QArguments, QException}
import collection.mutable.ListBuffer
import collection.JavaConversions._

/**
 * A basic job grapher.
 * Limitiations:
 *  - Only walks along graphs with rules that have a single input and a single output.
 */
class TreeJobGrapher extends JobGrapher {
  private val modelGraph = new SimpleDirectedGraph[QFile, QCommand](classOf[QCommand])

  private val rules = new ListBuffer[QRule]
  private var sourceFiles = List.empty[QFile]
  private var targetFiles = List.empty[QFile]

  override protected def createJobGraph() = {
    createModelGraph()
    
    var missingPaths = List.empty[Tuple2[QFile,QFile]]
    for (sourceFile <- sourceFiles) {
      for (targetFile <- targetFiles) {
        var shortestPath = BellmanFordShortestPath.findPathBetween(modelGraph, sourceFile, targetFile)
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
    for (rule <- rules) {
      if (rule.inputs.size != 1 || (rule.outputs.size != 1))
        throw new QException(this.getClass.getName + " can only process rules with a single input and a single output.  " +
                "inputs: " + rule.inputs + ", outputs: " + rule.outputs + ", command: " + rule.command)
      var source = rule.inputs.head
      var target = rule.outputs.head
      modelGraph.addVertex(source)
      modelGraph.addVertex(target)
      modelGraph.addEdge(source, target, rule.command)
    }
  }

  private def addPaths(shortestPath: java.util.List[QCommand], qArgs: QArguments) {
    for (inputFile <- qArgs.inputPaths.map(_.getCanonicalPath))
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
      val resourceSource = new MapResourceNode(Map(source.extension -> sourceFile))
      val resourceTarget = new MapResourceNode(Map(target.extension -> targetFile))
      val resourceEdge = new ExecEdge(command)
      jobGraph.addVertex(resourceSource)
      jobGraph.addVertex(resourceTarget)
      jobGraph.addEdge(resourceSource, resourceTarget, resourceEdge)
      sourceFile = targetFile
    }
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
