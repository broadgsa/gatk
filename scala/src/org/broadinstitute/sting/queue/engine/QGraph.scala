package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.function.{MappingFunction, CommandLineFunction, QFunction}
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.QException
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import java.io.File

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  var bsubWaitJobs = false
  val jobGraph = newGraph
  def numJobs = JavaConversions.asSet(jobGraph.edgeSet).filter(_.isInstanceOf[CommandLineFunction]).size

  /**
   * Adds a QScript created CommandLineFunction to the graph.
   * @param command Function to add to the graph.
   */
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

  /**
   * Checks the functions for missing values and the graph for cyclic dependencies and then runs the functions in the graph.
   */
  def run = {
    var isReady = true
    var totalMissingValues = 0
    for (function <- JavaConversions.asSet(jobGraph.edgeSet)) {
      function match {
        case cmd: CommandLineFunction =>
          val missingFieldValues = cmd.missingFields
          if (missingFieldValues.size > 0) {
            totalMissingValues += missingFieldValues.size
            logger.error("Missing %s values for function: %s".format(missingFieldValues.size, cmd.commandLine))
            for (missing <- missingFieldValues)
              logger.error("  " + missing)
          }
        case _ =>
      }
    }

    if (totalMissingValues > 0) {
      isReady = false
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

    if (totalMissingValues > 0) {
      logger.error("Total missing values: " + totalMissingValues)
    }

    if (isReady && this.dryRun) {
      logger.info("Dry run completed successfully!")
      logger.info("Re-run with \"-run\" to execute the functions.")
    }
  }

  /**
   * Creates a new graph where if new edges are needed (for cyclic dependency checking) they can be automatically created using a generic MappingFunction.
   * @return A new graph
   */
  private def newGraph = new SimpleDirectedGraph[QNode, QFunction](new EdgeFactory[QNode, QFunction] {
    def createEdge(input: QNode, output: QNode) = new MappingFunction(input.files, output.files)})

  /**
   * Adds a generic QFunction to the graph.
   * If the function is scatterable and the jobs request bsub, splits the job into parts and adds the parts instead.
   * @param f Generic QFunction to add to the graph.
   */
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

  /**
   * Checks to see if the set of files has more than one file and if so adds input mappings between the set and the individual files.
   * @param files Set to check.
   */
  private def addCollectionInputs(files: Set[File]): Unit = {
    if (files.size > 1)
      for (file <- files)
        addMappingEdge(Set(file), files)
  }

  /**
   * Checks to see if the set of files has more than one file and if so adds output mappings between the individual files and the set.
   * @param files Set to check.
   */
  private def addCollectionOutputs(files: Set[File]): Unit = {
    if (files.size > 1)
      for (file <- files)
        addMappingEdge(files, Set(file))
  }

  /**
   * Adds a directed graph edge between the input set and the output set if there isn't a direct relationship between the two nodes already.
   * @param input Input set of files.
   * @param output Output set of files.
   */
  private def addMappingEdge(input: Set[File], output: Set[File]) = {
    val hasEdge = input == output ||
            jobGraph.getEdge(QNode(input), QNode(output)) != null ||
            jobGraph.getEdge(QNode(output), QNode(input)) != null
    if (!hasEdge)
      addFunction(new MappingFunction(input, output))
  }

  /**
   * Returns true if the edge is an internal mapping edge.
   * @param edge Edge to check.
   * @return true if the edge is an internal mapping edge.
   */
  private def isMappingEdge(edge: QFunction) =
    edge.isInstanceOf[MappingFunction]

  /**
   * Returns true if the edge is mapping edge that is not needed because it does
   * not direct input or output from a user generated CommandLineFunction.
   * @param edge Edge to check.
   * @return true if the edge is not needed in the graph.
   */
  private def isFiller(edge: QFunction) = {
    if (isMappingEdge(edge)) {
      if (jobGraph.outgoingEdgesOf(jobGraph.getEdgeTarget(edge)).size == 0)
        true
      else if (jobGraph.incomingEdgesOf(jobGraph.getEdgeSource(edge)).size == 0)
        true
      else false
    } else false
  }

  /**
   * Returns true if the node is not connected to any edges.
   * @param node Node (set of files) to check
   * @return true if this set of files is not needed in the graph.
   */
  private def isOrphan(node: QNode) =
    (jobGraph.incomingEdgesOf(node).size + jobGraph.outgoingEdgesOf(node).size) == 0

  /**
   * Outputs the graph to a .dot file.
   * http://en.wikipedia.org/wiki/DOT_language
   * @param file Path to output the .dot file.
   */
  def renderToDot(file: java.io.File) = {
    val out = new java.io.FileWriter(file)

    // todo -- we need a nice way to visualize the key pieces of information about commands.  Perhaps a
    // todo -- visualizeString() command, or something that shows inputs / outputs
    val ve = new org.jgrapht.ext.EdgeNameProvider[QFunction] {
        def getEdgeName( function: QFunction ) = function.dotString
    }

    //val iterator = new TopologicalOrderIterator(qGraph.jobGraph)
    (new DOTExporter(new org.jgrapht.ext.IntegerNameProvider[QNode](), null, ve)).export(out, jobGraph)

    out.close
  }
}
