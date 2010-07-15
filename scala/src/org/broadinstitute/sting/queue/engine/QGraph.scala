package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.function.{MappingFunction, CommandLineFunction, QFunction}
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.queue.util.{CollectionUtils, Logging}
import org.broadinstitute.sting.queue.QException
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import org.broadinstitute.sting.queue.function.DispatchFunction
import org.broadinstitute.sting.queue.function.gatk.GatkFunction

class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  var bsubWaitJobs = false
  var properties = Map.empty[String, String]
  val jobGraph = newGraph
  def numJobs = JavaConversions.asSet(jobGraph.edgeSet).filter(_.isInstanceOf[CommandLineFunction]).size

  def add(command: CommandLineFunction) {
    addFunction(command)
  }

  /**
   * Looks through functions with multiple inputs and outputs and adds mapping functions for single inputs and outputs.
   */
  def fillIn = {
    // clone since edgeSet is backed by the graph
    for (function <- JavaConversions.asSet(jobGraph.edgeSet).clone) {
      addCollectionOutputs(function.outputs)
      addCollectionInputs(function.inputs)
    }

    var pruning = true
    while (pruning) {
      pruning = false
      val filler = jobGraph.edgeSet.filter(isFiller(_))
      if (filler.size > 0) {
        jobGraph.removeAllEdges(filler)
        pruning = true
      }
    }

    jobGraph.removeAllVertices(jobGraph.vertexSet.filter(isOrphan(_)))
  }

  def run = {
    var isReady = true
    for (function <- JavaConversions.asSet(jobGraph.edgeSet)) {
      function match {
        case cmd: CommandLineFunction =>
          val missingValues = cmd.missingValues
          if (missingValues.size > 0) {
            isReady = false
            logger.error("Missing values for function: %s".format(cmd.commandLine))
            for (missing <- missingValues)
              logger.error("  " + missing)
          }
        case _ =>
      }
    }

    val detector = new CycleDetector(jobGraph)
    if (detector.detectCycles) {
      logger.error("Cycles were detected in the graph:")
      for (cycle <- detector.findCycles)
        logger.error("  " + cycle)
      isReady = false
    }

    if (isReady || this.dryRun)
      (new TopologicalJobScheduler(this) with LsfJobRunner).runJobs
  }

  private def newGraph = new SimpleDirectedGraph[QNode, QFunction](new EdgeFactory[QNode, QFunction] {
    def createEdge(input: QNode, output: QNode) = new MappingFunction(input.items, output.items)})

  private def addFunction(f: QFunction): Unit = {
    try {
      f.freeze

      f match {
        case scatterGather: ScatterGatherableFunction if (bsubAllJobs && scatterGather.scatterGatherable) =>
          val functions = scatterGather.generateFunctions()
          if (logger.isTraceEnabled)
            logger.trace("Scattered into %d parts: %s".format(functions.size, functions))
          functions.foreach(addFunction(_))
        case _ =>
          val inputs = QNode(f.inputs)
          val outputs = QNode(f.outputs)
          val newSource = jobGraph.addVertex(inputs)
          val newTarget = jobGraph.addVertex(outputs)
          val removedEdges = jobGraph.removeAllEdges(inputs, outputs)
          val added = jobGraph.addEdge(inputs, outputs, f)
          if (logger.isTraceEnabled) {
            logger.trace("Mapped from:   " + inputs)
            logger.trace("Mapped to:     " + outputs)
            logger.trace("Mapped via:    " + f)
            logger.trace("Removed edges: " + removedEdges)
            logger.trace("New source?:   " + newSource)
            logger.trace("New target?:   " + newTarget)
            logger.trace("")
          }
      }
    } catch {
      case e: Exception =>
        throw new QException("Error adding function: " + f, e)
    }
  }

  private def addCollectionInputs(value: Any): Unit = {
    CollectionUtils.foreach(value, (item, collection) =>
      addMappingEdge(item, collection))
  }

  private def addCollectionOutputs(value: Any): Unit = {
    CollectionUtils.foreach(value, (item, collection) =>
      addMappingEdge(collection, item))
  }

  private def addMappingEdge(input: Any, output: Any) = {
    val inputSet = asSet(input)
    val outputSet = asSet(output)
    val hasEdge = inputSet == outputSet ||
            jobGraph.getEdge(QNode(inputSet), QNode(outputSet)) != null ||
            jobGraph.getEdge(QNode(outputSet), QNode(inputSet)) != null
    if (!hasEdge)
      addFunction(new MappingFunction(inputSet, outputSet))
  }

  private def asSet(value: Any): Set[Any] = if (value.isInstanceOf[Set[_]]) value.asInstanceOf[Set[Any]] else Set(value)

  private def isMappingEdge(edge: QFunction) =
    edge.isInstanceOf[MappingFunction]

  private def isFiller(edge: QFunction) = {
    if (isMappingEdge(edge)) {
      if (jobGraph.outgoingEdgesOf(jobGraph.getEdgeTarget(edge)).size == 0)
        true
      else if (jobGraph.incomingEdgesOf(jobGraph.getEdgeSource(edge)).size == 0)
        true
      else false
    } else false
  }

  private def isOrphan(node: QNode) =
    (jobGraph.incomingEdgesOf(node).size + jobGraph.outgoingEdgesOf(node).size) == 0

  def renderToDot(file: java.io.File) = {
    val out = new java.io.FileWriter(file)

    // todo -- we need a nice way to visualize the key pieces of information about commands.  Perhaps a
    // todo -- visualizeString() command, or something that shows inputs / outputs
    val ve = new org.jgrapht.ext.EdgeNameProvider[QFunction] {
        def getEdgeName( function: QFunction ) = function match {
            case f: DispatchFunction => f.jobName + " => " + f.commandLine
            case _ => ""
        }
    }

    //val iterator = new TopologicalOrderIterator(qGraph.jobGraph)
    (new DOTExporter(new org.jgrapht.ext.IntegerNameProvider[QNode](), null, ve)).export(out, jobGraph)

    out.close
  }
}
