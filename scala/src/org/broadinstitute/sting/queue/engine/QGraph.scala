package org.broadinstitute.sting.queue.engine

import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.JavaConversions._
import scala.collection.immutable.ListMap
import org.broadinstitute.sting.queue.function.{MappingFunction, CommandLineFunction, QFunction}
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.queue.util.{CollectionUtils, Logging}
import org.broadinstitute.sting.queue.QException
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory

class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  val jobGraph = newGraph
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

      for ((name, input) <- inputs) {
        addCollectionInputs(name, input)
        if (inputs.size > 1)
          addMappingEdge(ListMap(name -> input), inputs)
      }

      for ((name, output) <- outputs) {
        addCollectionOutputs(name, output)
        if (outputs.size > 1)
          addMappingEdge(outputs, ListMap(name -> output))
      }
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
    def createEdge(input: QNode, output: QNode) = new MappingFunction(input.valueMap, output.valueMap)})

  private def add(f: QFunction, replace: Boolean): Unit = {
    try {
      f.freeze

      f match {
        case scatterGather: ScatterGatherableFunction if (bsubAllJobs && scatterGather.scatterGatherable) =>
          val functions = scatterGather.generateFunctions()
          if (logger.isTraceEnabled)
            logger.trace("Scattered into %d parts: %s".format(functions.size, functions))
          functions.foreach(add(_))
        case _ =>
          val inputs = QNode(f.inputs.values.filter(_ != null).toSet)
          val outputs = QNode(f.outputs.values.filter(_ != null).toSet)
          val newSource = jobGraph.addVertex(inputs)
          val newTarget = jobGraph.addVertex(outputs)
          val removedEdges = if (replace) jobGraph.removeAllEdges(inputs, outputs) else Nil
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

  private def addCollectionInputs(name: String, value: Any): Unit = {
    CollectionUtils.foreach(value, (item, collection) =>
      addMappingEdge(ListMap(name -> item), ListMap(name -> collection)))
  }

  private def addCollectionOutputs(name: String, value: Any): Unit = {
    CollectionUtils.foreach(value, (item, collection) =>
      addMappingEdge(ListMap(name -> collection), ListMap(name -> item)))
  }

  private def addMappingEdge(input: ListMap[String, Any], output: ListMap[String, Any]) =
    add(new MappingFunction(input, output), false)

  private def isMappingEdge(edge: QFunction) =
    edge.isInstanceOf[MappingFunction]

  private def isFiller(edge: QFunction) = {
    if (isMappingEdge(edge)) {
      val source = jobGraph.getEdgeSource(edge)
      val target = jobGraph.getEdgeTarget(edge)
      if (jobGraph.outgoingEdgesOf(target).size == 0 || jobGraph.incomingEdgesOf(source).size == 0)
        true
      else if (isLoopback(source) || isLoopback(target))
        true
      else false
    } else false
  }

  private def isLoopback(node: QNode) = {
    var loopback = false
    val incoming = jobGraph.incomingEdgesOf(node)
    val outgoing = jobGraph.outgoingEdgesOf(node)
    if (incoming.size == 1 && outgoing.size == 1)
      if (isMappingEdge(incoming.head) && isMappingEdge(outgoing.head))
        if (jobGraph.getEdgeSource(incoming.head) == jobGraph.getEdgeTarget(outgoing.head))
          loopback = true
    loopback
  }

  private def isOrphan(node: QNode) =
    (jobGraph.incomingEdgesOf(node).size + jobGraph.outgoingEdgesOf(node).size) == 0
}
