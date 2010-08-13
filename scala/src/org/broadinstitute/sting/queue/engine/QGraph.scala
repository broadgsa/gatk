package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.queue.util.Logging
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import java.io.File
import org.jgrapht.event.{TraversalListenerAdapter, EdgeTraversalEvent}
import org.broadinstitute.sting.queue.{QSettings, QException}
import org.broadinstitute.sting.queue.function.{DispatchWaitFunction, MappingFunction, CommandLineFunction, QFunction}

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  var bsubWaitJobs = false
  var skipUpToDateJobs = false
  var dotFile: File = _
  var expandedDotFile: File = _
  var qSettings: QSettings = _
  var debugMode = false
  private val jobGraph = newGraph

  /**
   * Adds a QScript created CommandLineFunction to the graph.
   * @param command Function to add to the graph.
   */
  def add(command: CommandLineFunction) {
    addFunction(command)
  }

  /**
   * Checks the functions for missing values and the graph for cyclic dependencies and then runs the functions in the graph.
   */
  def run = {
    fill
    if (dotFile != null)
      renderToDot(dotFile)
    var numMissingValues = validate

    if (numMissingValues == 0 && bsubAllJobs) {
      logger.debug("Scatter gathering jobs.")
      var scatterGathers = List.empty[ScatterGatherableFunction]
      loop({
        case scatterGather: ScatterGatherableFunction if (scatterGather.scatterGatherable) =>
          scatterGathers :+= scatterGather
      })

      var addedFunctions = List.empty[CommandLineFunction]
      for (scatterGather <- scatterGathers) {
        val functions = scatterGather.generateFunctions()
        if (this.debugMode)
          logger.debug("Scattered into %d parts: %n%s".format(functions.size, functions.mkString("%n".format())))
        addedFunctions ++= functions
      }

      this.jobGraph.removeAllEdges(scatterGathers)
      prune
      addedFunctions.foreach(this.addFunction(_))

      fill
      val scatterGatherDotFile = if (expandedDotFile != null) expandedDotFile else dotFile
      if (scatterGatherDotFile != null)
        renderToDot(scatterGatherDotFile)
      numMissingValues = validate
    }

    val isReady = numMissingValues == 0

    if (isReady || this.dryRun)
      runJobs

    if (numMissingValues > 0) {
      logger.error("Total missing values: " + numMissingValues)
    }

    if (isReady && this.dryRun) {
      logger.info("Dry run completed successfully!")
      logger.info("Re-run with \"-run\" to execute the functions.")
    }
  }

  /**
   * Walks up the graph looking for the previous LsfJobs.
   * @param function Function to examine for a previous command line job.
   * @param qGraph The graph that contains the jobs.
   * @return A list of prior jobs.
   */
  def previousJobs(function: QFunction) : List[CommandLineFunction] = {
    var previous = List.empty[CommandLineFunction]

    val source = this.jobGraph.getEdgeSource(function)
    for (incomingEdge <- this.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

      // Stop recursing when we find a job along the edge and return its job id
        case commandLineFunction: CommandLineFunction => previous :+= commandLineFunction

        // For any other type of edge find the LSF jobs preceding the edge
        case qFunction: QFunction => previous ++= previousJobs(qFunction)
      }
    }
    previous
  }

  /**
   * Fills in the graph using mapping functions, then removes out of date
   * jobs, then cleans up mapping functions and nodes that aren't need.
   */
  private def fill = {
    fillIn
    if (skipUpToDateJobs)
      removeUpToDate
    prune
  }

  /**
   * Looks through functions with multiple inputs and outputs and adds mapping functions for single inputs and outputs.
   */
  private def fillIn = {
    // clone since edgeSet is backed by the graph
    JavaConversions.asSet(jobGraph.edgeSet).clone.foreach {
      case cmd: CommandLineFunction => {
        addCollectionOutputs(cmd.outputs)
        addCollectionInputs(cmd.inputs)
      }
      case map: MappingFunction => /* do nothing for mapping functions */
    }
  }

  /**
   * Removes functions that are up to date.
   */
  private def removeUpToDate = {
    var upToDateJobs = Set.empty[CommandLineFunction]
    loop({
      case f if (upToDate(f, upToDateJobs)) => {
        logger.info("Skipping command because it is up to date: %n%s".format(f.commandLine))
        upToDateJobs += f
      }
    })
    for (upToDateJob <- upToDateJobs)
      jobGraph.removeEdge(upToDateJob)
  }

  /**
   * Returns true if the all previous functions in the graph are up to date, and the function is up to date.
   */
  private def upToDate(commandLineFunction: CommandLineFunction, upToDateJobs: Set[CommandLineFunction]) = {
    this.previousJobs(commandLineFunction).forall(upToDateJobs.contains(_)) && commandLineFunction.upToDate
  }

  /**
   *  Removes mapping edges that aren't being used, and nodes that don't belong to anything.
   */
  private def prune = {
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
   * Validates that the functions in the graph have no missing values and that there are no cycles.
   * @return Number of missing values.
   */
  private def validate = {
    var numMissingValues = 0
    JavaConversions.asSet(jobGraph.edgeSet).foreach {
      case cmd: CommandLineFunction =>
        val missingFieldValues = cmd.missingFields
        if (missingFieldValues.size > 0) {
          numMissingValues += missingFieldValues.size
          logger.error("Missing %s values for function: %s".format(missingFieldValues.size, cmd.commandLine))
          for (missing <- missingFieldValues)
            logger.error("  " + missing)
        }
      case map: MappingFunction => /* do nothing for mapping functions */
    }

    val detector = new CycleDetector(jobGraph)
    if (detector.detectCycles) {
      logger.error("Cycles were detected in the graph:")
      for (cycle <- detector.findCycles)
        logger.error("  " + cycle)
      throw new QException("Cycles were detected in the graph.")
    }

    numMissingValues
  }

  /**
   * Runs the jobs by traversing the graph.
   */
  private def runJobs = {
    val runner = if (bsubAllJobs) new LsfJobRunner else new ShellJobRunner

    val numJobs = JavaConversions.asSet(jobGraph.edgeSet).filter(_.isInstanceOf[CommandLineFunction]).size

    logger.info("Number of jobs: %s".format(numJobs))
    if (this.debugMode) {
      val numNodes = jobGraph.vertexSet.size
      logger.debug("Number of nodes: %s".format(numNodes))
    }
    var numNodes = 0

    loop(
      edgeFunction = { case f => runner.run(f, this) },
      nodeFunction = {
        case node => {
          if (this.debugMode)
            logger.debug("Visiting: " + node)
          numNodes += 1
        }
      })

    if (this.debugMode)
      logger.debug("Done walking %s nodes.".format(numNodes))

    if (bsubAllJobs && bsubWaitJobs) {
      logger.info("Waiting for jobs to complete.")
      val wait = new DispatchWaitFunction
      wait.qSettings = this.qSettings
      wait.freeze
      runner.run(wait, this)
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
   * @param f Generic QFunction to add to the graph.
   */
  private def addFunction(f: QFunction): Unit = {
    try {
      f match {
        case cmd: CommandLineFunction => cmd.qSettings = this.qSettings
        case map: MappingFunction => /* do nothing for mapping functions */
      }
      f.freeze
      val inputs = QNode(f.inputs)
      val outputs = QNode(f.outputs)
      val newSource = jobGraph.addVertex(inputs)
      val newTarget = jobGraph.addVertex(outputs)
      val removedEdges = jobGraph.removeAllEdges(inputs, outputs)
      val added = jobGraph.addEdge(inputs, outputs, f)
      if (this.debugMode) {
        logger.debug("Mapped from:   " + inputs)
        logger.debug("Mapped to:     " + outputs)
        logger.debug("Mapped via:    " + f)
        logger.debug("Removed edges: " + removedEdges)
        logger.debug("New source?:   " + newSource)
        logger.debug("New target?:   " + newTarget)
        logger.debug("")
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
   * Utility function for looping over the internal graph and running functions.
   * @param edgeFunction Optional function to run for each edge visited.
   * @param nodeFunction Optional function to run for each node visited.
   */
  private def loop(edgeFunction: PartialFunction[CommandLineFunction, Unit] = null, nodeFunction: PartialFunction[QNode, Unit] = null) = {
    val iterator = new TopologicalOrderIterator(this.jobGraph)
    iterator.addTraversalListener(new TraversalListenerAdapter[QNode, QFunction] {
      override def edgeTraversed(event: EdgeTraversalEvent[QNode, QFunction]) = event.getEdge match {
        case cmd: CommandLineFunction => if (edgeFunction != null && edgeFunction.isDefinedAt(cmd)) edgeFunction(cmd)
        case map: MappingFunction => /* do nothing for mapping functions */
      }
    })
    iterator.foreach(node => if (nodeFunction != null && nodeFunction.isDefinedAt(node)) nodeFunction(node))
  }

  /**
   * Outputs the graph to a .dot file.
   * http://en.wikipedia.org/wiki/DOT_language
   * @param file Path to output the .dot file.
   */
  private def renderToDot(file: java.io.File) = {
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
