package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions._
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import java.io.File
import org.jgrapht.event.{TraversalListenerAdapter, EdgeTraversalEvent}
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.queue.function.{InProcessFunction, CommandLineFunction, QFunction}
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, GatherFunction, ScatterGatherableFunction}
import org.apache.commons.lang.StringUtils
import org.broadinstitute.sting.queue.util._

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var settings: QGraphSettings = _
  var debugMode = false

  private def dryRun = !settings.run
  private var numMissingValues = 0
  private val jobGraph = newGraph
  private var running = true
  private val runningLock = new Object
  private val nl = "%n".format()

  private val inProcessManager = new InProcessJobManager
  private var commandLineManager: JobManager[CommandLineFunction, _<:JobRunner[CommandLineFunction]] = _

  /**
   * Adds a QScript created CommandLineFunction to the graph.
   * @param command Function to add to the graph.
   */
  def add(command: QFunction) {
    try {
      runningLock.synchronized {
        if (running) {
          command.qSettings = settings.qSettings
          command.freeze
          val inputs = getQNode(command.inputs.toList.sortWith(_.compareTo(_) < 0))
          val outputs = getQNode(command.outputs.toList.sortWith(_.compareTo(_) < 0))
          addEdge(new FunctionEdge(command, inputs, outputs))
        }
      }
    } catch {
      case e: Exception =>
        throw new QException("Error adding function: " + command, e)
    }
  }

  /**
   * Checks the functions for missing values and the graph for cyclic dependencies and then runs the functions in the graph.
   */
  def run() {
    runningLock.synchronized {
      if (running) {
        IOUtils.checkTempDir(settings.qSettings.tempDirectory)
        fillGraph
        val isReady = numMissingValues == 0

        if (this.jobGraph.edgeSet.isEmpty) {
          logger.warn("Nothing to run! Were any Functions added?");
        } else if (settings.getStatus) {
          logger.info("Checking pipeline status.")
          logStatus()
        } else if (this.dryRun) {
          dryRunJobs()
        } else if (isReady) {
          logger.info("Running jobs.")
          runJobs()
        }

        if (numMissingValues > 0) {
          logger.error("Total missing values: " + numMissingValues)
        }

        if (running && isReady && this.dryRun) {
          logger.info("Dry run completed successfully!")
          logger.info("Re-run with \"-run\" to execute the functions.")
        }
      }
    }
  }

  private def fillGraph {
    logger.info("Generating graph.")
    fill
    if (settings.dotFile != null)
      renderToDot(settings.dotFile)
    validate()

    if (running && numMissingValues == 0 && settings.bsubAllJobs) {
      logger.info("Generating scatter gather jobs.")
      val scatterGathers = jobGraph.edgeSet.filter(edge => scatterGatherable(edge))

      var addedFunctions = List.empty[QFunction]
      for (scatterGather <- scatterGathers) {
        val functions = scatterGather.asInstanceOf[FunctionEdge]
                .function.asInstanceOf[ScatterGatherableFunction]
                .generateFunctions()
        if (this.debugMode)
          logger.debug("Scattered into %d parts: %n%s".format(functions.size, functions.mkString(nl)))
        addedFunctions ++= functions
      }

      logger.info("Removing original jobs.")
      this.jobGraph.removeAllEdges(scatterGathers)
      prune()

      logger.info("Adding scatter gather jobs.")
      addedFunctions.foreach(function => if (running) this.add(function))

      logger.info("Regenerating graph.")
      fill
      val scatterGatherDotFile = if (settings.expandedDotFile != null) settings.expandedDotFile else settings.dotFile
      if (scatterGatherDotFile != null)
        renderToDot(scatterGatherDotFile)
      validate()
    }
  }

  private def scatterGatherable(edge: QEdge) = {
    edge match {
      case functionEdge: FunctionEdge => {
        functionEdge.function match {
          case scatterGather: ScatterGatherableFunction if (scatterGather.scatterGatherable) => true
          case _ => false
        }
      }
      case _ => false
    }
  }

  /**
   * Walks up the graph looking for the previous command line edges.
   * @param function Function to examine for a previous command line job.
   * @param qGraph The graph that contains the jobs.
   * @return A list of prior jobs.
   */
  private def previousFunctions(edge: QEdge): List[FunctionEdge] = {
    var previous = List.empty[FunctionEdge]

    val source = this.jobGraph.getEdgeSource(edge)
    for (incomingEdge <- this.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

        // Stop recursing when we find a job along the edge and return its job id
        case functionEdge: FunctionEdge => previous :+= functionEdge

        // For any other type of edge find the jobs preceding the edge
        case edge: QEdge => previous ++= previousFunctions(edge)
      }
    }
    previous
  }

  /**
   * Fills in the graph using mapping functions, then removes out of date
   * jobs, then cleans up mapping functions and nodes that aren't need.
   */
  private def fill() {
    fillIn()
    prune()
  }

  /**
   * Looks through functions with multiple inputs and outputs and adds mapping functions for single inputs and outputs.
   */
  private def fillIn() {
    // clone since edgeSet is backed by the graph
    asScalaSet(jobGraph.edgeSet).clone.foreach(edge => {
      if (running) edge match {
        case cmd: FunctionEdge => {
          addCollectionOutputs(cmd.outputs)
          addCollectionInputs(cmd.inputs)
        }
        case map: MappingEdge => /* do nothing for mapping edges */
      }
    })
  }

  private def getReadyJobs = {
    jobGraph.edgeSet.filter{
      case f: FunctionEdge =>
        this.previousFunctions(f).forall(_.status == RunnerStatus.DONE) && f.status == RunnerStatus.PENDING
      case _ => false
    }.toList.asInstanceOf[List[FunctionEdge]].sortWith(compare(_,_))
  }

  private def getRunningJobs = {
    jobGraph.edgeSet.filter{
      case f: FunctionEdge => f.status == RunnerStatus.RUNNING
      case _ => false
    }.toList.asInstanceOf[List[FunctionEdge]].sortWith(compare(_,_))
  }

  /**
   *  Removes mapping edges that aren't being used, and nodes that don't belong to anything.
   */
  private def prune() {
    var pruning = true
    while (pruning) {
      pruning = false
      val filler = jobGraph.edgeSet.filter(isFiller(_))
      if (filler.size > 0) {
        jobGraph.removeAllEdges(filler)
        pruning = running
      }
    }

    if (running)
      jobGraph.removeAllVertices(jobGraph.vertexSet.filter(isOrphan(_)))
  }

  /**
   * Validates that the functions in the graph have no missing values and that there are no cycles.
   */
  private def validate() {
    asScalaSet(jobGraph.edgeSet).foreach(
      edge =>
        if (running) edge match
        {
          case cmd: FunctionEdge =>
            val missingFieldValues = cmd.function.missingFields
            if (missingFieldValues.size > 0) {
              numMissingValues += missingFieldValues.size
              logger.error("Missing %s values for function: %s".format(missingFieldValues.size, cmd.function.description))
              for (missing <- missingFieldValues)
                logger.error("  " + missing)
            }
          case map: MappingEdge => /* do nothing for mapping edges */
        }
    )

    val detector = new CycleDetector(jobGraph)
    if (detector.detectCycles) {
      logger.error("Cycles were detected in the graph:")
      for (cycle <- detector.findCycles)
        logger.error("  " + cycle)
      throw new QException("Cycles were detected in the graph.")
    }
  }

  /**
   * Dry-runs the jobs by traversing the graph.
   */
  private def dryRunJobs() {
    updateGraphStatus(false)
    var readyJobs = getReadyJobs
    while (running && readyJobs.size > 0) {
      logger.debug("+++++++")
      readyJobs.foreach(edge => {
        if (running) {
          logEdge(edge)
          edge.markAsDone
        }
      })
      readyJobs = getReadyJobs
    }
  }

  private def logEdge(edge: FunctionEdge) {
    logger.info("-------")
    if (logger.isDebugEnabled) {
      logger.debug("Inputs: " + edge.inputs)
    }
    logger.info(StringUtils.capitalize(edge.status.toString) + ": " + edge.function.description)
    if (logger.isDebugEnabled)
      logger.debug(edge.function.commandDirectory + " > " + edge.function.description)
    logger.info("Log: " + edge.function.jobOutputFile.getAbsolutePath)
    if (edge.function.jobErrorFile != null)
      logger.info("Error: " + edge.function.jobErrorFile.getAbsolutePath)
  }

  /**
   * Logs job statuses by traversing the graph and looking for status-related files
   */
  private def logStatus() {
    updateGraphStatus(false)
    doStatus(status => logger.info(status))
  }

  /**
   * Runs the jobs by traversing the graph.
   */
  private def runJobs() {
    try {
      if (settings.bsubAllJobs)
        commandLineManager = new Lsf706JobManager
      else
        commandLineManager = new ShellJobManager

      if (settings.startFromScratch) {
        logger.info("Removing outputs from previous runs.")
        foreachFunction(_.resetToPending(true))
      } else
        updateGraphStatus(true)

      var readyJobs = getReadyJobs
      var runningJobs = Set.empty[FunctionEdge]
      while (running && readyJobs.size + runningJobs.size > 0) {
        var exitedJobs = List.empty[FunctionEdge]
        var failedJobs = List.empty[FunctionEdge]

        runningJobs.foreach(runner => runner.status match {
          case RunnerStatus.RUNNING => /* do nothing while still running */
          case RunnerStatus.FAILED => exitedJobs :+= runner; failedJobs :+= runner
          case RunnerStatus.DONE => exitedJobs :+= runner
        })
        exitedJobs.foreach(runner => runningJobs -= runner)

        readyJobs.foreach(f => {
          if (running) {
            f.runner = newRunner(f.function)
            f.start()
            f.status match {
              case RunnerStatus.RUNNING => runningJobs += f
              case RunnerStatus.FAILED => failedJobs :+= f
              case RunnerStatus.DONE => /* do nothing and move on */
            }
          }
        })

        if (failedJobs.size > 0) {
          emailFailedJobs(failedJobs)
          checkRetryJobs(failedJobs)
        }

        if (readyJobs.size == 0 && runningJobs.size > 0)
          Thread.sleep(30000L)
        readyJobs = getReadyJobs
      }

      deleteIntermediateOutputs()
    } catch {
      case e =>
        logger.error("Uncaught error running jobs.", e)
        throw e
    } finally {
      emailStatus()
    }
  }

  /**
   * Updates the status of edges in the graph.
   * @param cleanOutputs If true will delete outputs when setting edges to pending.
   */
  private def updateGraphStatus(cleanOutputs: Boolean) {
    traverseFunctions(edge => checkDone(edge, cleanOutputs))
  }

  /**
   * Checks if an edge is done or if it's an intermediate edge if it can be skipped.
   * This function may modify previous edges if it discovers that the edge passed in
   * is dependent jobs that were previously marked as skipped.
   * @param edge Edge to check to see if it's done or can be skipped.
   * @param cleanOutputs If true will delete outputs when setting edges to pending.
   */
  private def checkDone(edge: FunctionEdge, cleanOutputs: Boolean) {
    if (edge.function.isIntermediate) {
      // By default we do not need to run intermediate edges.
      // Mark any intermediate edges as skipped, if they're not already done.
      if (edge.status != RunnerStatus.DONE)
        edge.markAsSkipped()
    } else {
      val previous = this.previousFunctions(edge)
      val isDone = edge.status == RunnerStatus.DONE &&
              previous.forall(edge => edge.status == RunnerStatus.DONE || edge.status == RunnerStatus.SKIPPED)
      if (!isDone) {
        edge.resetToPending(cleanOutputs)
        resetPreviousSkipped(edge, previous, cleanOutputs)
      }
    }
  }

  /**
   * From the previous edges, resets any that are marked as skipped to pending.
   * If those that are reset have skipped edges, those skipped edges are recursively also set
   * to pending.
   * @param edge Dependent edge.
   * @param previous Previous edges that provide inputs to edge.
   * @param cleanOutputs If true will clean up the output files when resetting skipped jobs to pending.
   */
  private def resetPreviousSkipped(edge: FunctionEdge, previous: List[FunctionEdge], cleanOutputs: Boolean) {
    for (previousEdge <- previous.filter(_.status == RunnerStatus.SKIPPED)) {
      previousEdge.resetToPending(cleanOutputs)
      resetPreviousSkipped(previousEdge, this.previousFunctions(previousEdge), cleanOutputs)
    }
  }

  private def newRunner(f: QFunction) = {
    f match {
      case cmd: CommandLineFunction =>
        commandLineManager.create(cmd)
      case inProc: InProcessFunction =>
        inProcessManager.create(inProc)
      case _ =>
        throw new QException("Unexpected function: " + f)
    }
  }

  private def emailFailedJobs(failed: List[FunctionEdge]) {
    if (running && settings.statusEmailTo.size > 0) {
      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.subject = "Queue function: Failure: " + settings.qSettings.jobNamePrefix
      addFailedFunctions(emailMessage, failed)
      emailMessage.trySend(settings.qSettings.emailSettings)
    }
  }

  private def checkRetryJobs(failed: List[FunctionEdge]) {
    if (settings.retries > 0) {
      for (failedJob <- failed) {
        if (failedJob.function.jobRestartable && failedJob.retries < settings.retries) {
          failedJob.retries += 1
          failedJob.resetToPending(true)
          logger.info("Reset for retry attempt %d of %d: %s".format(
            failedJob.retries, settings.retries, failedJob.function.description))
        } else {
          logger.info("Giving up after retrying %d times: %s".format(
            settings.retries, failedJob.function.description))
        }
      }
    }
  }

  private def emailStatus() {
    if (running && settings.statusEmailTo.size > 0) {
      var failed = List.empty[FunctionEdge]
      foreachFunction(edge => {
        if (edge.status == RunnerStatus.FAILED) {
          failed :+= edge
        }
      })

      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.body = getStatus + nl
      if (failed.size == 0) {
        emailMessage.subject = "Queue run: Success: " + settings.qSettings.jobNamePrefix
      } else {
        emailMessage.subject = "Queue run: Failure: " + settings.qSettings.jobNamePrefix
        addFailedFunctions(emailMessage, failed)
      }
      emailMessage.trySend(settings.qSettings.emailSettings)
    }
  }

  private def addFailedFunctions(emailMessage: EmailMessage, failed: List[FunctionEdge]) {
    val logs = failed.flatMap(edge => logFiles(edge))

    if (emailMessage.body == null)
      emailMessage.body = ""
    emailMessage.body += """
    |Failed functions:
    |
    |%s
    |
    |Logs:
    |%s%n
    |""".stripMargin.trim.format(
      failed.map(edge => failedDescription(edge)).mkString(nl+nl),
      logs.map(_.getAbsolutePath).mkString(nl))

    emailMessage.attachments = logs
  }

  private def failedDescription(failed: FunctionEdge) = {
    var description = new StringBuilder
    if (settings.retries > 0)
      description.append("Attempt %d of %d.%n".format(failed.retries + 1, settings.retries + 1))
    description.append(failed.function.description)
    description.toString
  }

  private def logFiles(edge: FunctionEdge) = {
    var failedOutputs = List.empty[File]
    failedOutputs :+= edge.function.jobOutputFile
    if (edge.function.jobErrorFile != null)
      failedOutputs :+= edge.function.jobErrorFile
    failedOutputs.filter(file => file != null && file.exists)
  }

  /**
   * Tracks analysis status.
   */
  private class AnalysisStatus(val analysisName: String) {
    var status = RunnerStatus.PENDING
    var scatter = new ScatterGatherStatus
    var gather = new ScatterGatherStatus
  }

  /**
   * Tracks scatter gather status.
   */
  private class ScatterGatherStatus {
    var total = 0
    var done = 0
    var failed = 0
    var skipped = 0
  }

  /**
   * Gets job statuses by traversing the graph and looking for status-related files
   */
  private def getStatus = {
    val buffer = new StringBuilder
    doStatus(status => buffer.append(status).append(nl))
    buffer.toString
  }

  /**
   * Gets job statuses by traversing the graph and looking for status-related files
   */
  private def doStatus(statusFunc: String => Unit) = {
    var statuses = List.empty[AnalysisStatus]
    var maxWidth = 0
    foreachFunction(edge => {
      val name = edge.function.analysisName
      if (name != null) {
        updateAnalysisStatus(statuses.find(_.analysisName == name) match {
          case Some(status) => status
          case None =>
            val status = new AnalysisStatus(name)
            maxWidth = maxWidth max name.length
            statuses :+= status
            status
        }, edge)
      }
    })

    statuses.foreach(status => {
      val sgTotal = status.scatter.total + status.gather.total
      val sgDone = status.scatter.done + status.gather.done
      val sgFailed = status.scatter.failed + status.gather.failed
      val sgSkipped = status.scatter.skipped + status.gather.skipped
      val gatherTotal = status.gather.total
      val gatherDone = status.gather.done
      if (sgTotal > 0) {
        var sgStatus = RunnerStatus.PENDING
        if (sgFailed > 0)
          sgStatus = RunnerStatus.FAILED
        else if (gatherDone == gatherTotal)
          sgStatus = RunnerStatus.DONE
        else if (sgDone + sgSkipped == sgTotal)
          sgStatus = RunnerStatus.SKIPPED
        else if (sgDone > 0)
          sgStatus = RunnerStatus.RUNNING
        status.status = sgStatus
      }

      var info = ("%-" + maxWidth + "s [%s]")
              .format(status.analysisName, StringUtils.center(status.status.toString, 7))
      if (status.scatter.total + status.gather.total > 1) {
        info += formatSGStatus(status.scatter, "s")
        info += formatSGStatus(status.gather, "g")
      }
      statusFunc(info)
    })
  }

  /**
   * Updates a status map with scatter/gather status information (e.g. counts)
   */
  private def updateAnalysisStatus(stats: AnalysisStatus, edge: FunctionEdge) {
    if (edge.function.isInstanceOf[GatherFunction]) {
      updateSGStatus(stats.gather, edge)
    } else if (edge.function.isInstanceOf[CloneFunction]) {
      updateSGStatus(stats.scatter, edge)
    } else {
      stats.status = edge.status
    }
  }

  private def updateSGStatus(stats: ScatterGatherStatus, edge: FunctionEdge) {
    stats.total += 1
    edge.status match {
      case RunnerStatus.DONE => stats.done += 1
      case RunnerStatus.FAILED => stats.failed += 1
      case RunnerStatus.SKIPPED => stats.skipped += 1
      /* can't tell the difference between pending and running right now! */
      case RunnerStatus.PENDING =>
      case RunnerStatus.RUNNING =>
    }
  }

  /**
   * Formats a status into nice strings
   */
  private def formatSGStatus(stats: ScatterGatherStatus, prefix: String) = {
    " %s:%dt/%dd/%df".format(
      prefix, stats.total, stats.done, stats.failed)
  }

  /**
   *   Creates a new graph where if new edges are needed (for cyclic dependency checking) they can be automatically created using a generic MappingFunction.
   * @return A new graph
   */
  private def newGraph = new SimpleDirectedGraph[QNode, QEdge](new EdgeFactory[QNode, QEdge] {
    def createEdge(input: QNode, output: QNode) = new MappingEdge(input, output)})

  private var nextNodeId = 0
  private def getQNode(files: List[File]) = {
    jobGraph.vertexSet.find(node => node.files == files) match {
      case Some(node) =>
        node
      case None =>
        if (nextNodeId % 100 == 0)
          logger.debug("adding QNode: " + nextNodeId)
        val node = new QNode(nextNodeId, files)
        nextNodeId += 1
        jobGraph.addVertex(node)
        node
    }
  }

  private def addEdge(edge: QEdge) {
    jobGraph.removeAllEdges(edge.inputs, edge.outputs)
    jobGraph.addEdge(edge.inputs, edge.outputs, edge)
  }

  /**
   * Adds input mappings between the node's files and the individual files.
   * @param inputs Input node.
   */
  private def addCollectionInputs(inputs: QNode) {
    if (inputs.files.size > 1)
      for (file <- inputs.files) {
        if (running) {
          val input = getQNode(List(file))
          if (!jobGraph.containsEdge(input, inputs))
            addEdge(new MappingEdge(input, inputs))
        }
      }
  }

  /**
   * Adds output mappings between the node's files and the individual files.
   * @param outputs Output node.
   */
  private def addCollectionOutputs(outputs: QNode) {
    if (outputs.files.size > 1)
      for (file <- outputs.files) {
        if (running) {
          val output = getQNode(List(file))
          if (!jobGraph.containsEdge(outputs, output))
            addEdge(new MappingEdge(outputs, output))
        }
      }
  }

  /**
   * Returns true if the edge is mapping edge that is not needed because it does
   * not direct input or output from a user generated CommandLineFunction.
   * @param edge Edge to check.
   * @return true if the edge is not needed in the graph.
   */
  private def isFiller(edge: QEdge) = {
    edge match {
      case mapping: MappingEdge =>
        jobGraph.outgoingEdgesOf(jobGraph.getEdgeTarget(edge)).size == 0 &&
          jobGraph.incomingEdgesOf(jobGraph.getEdgeSource(edge)).size == 0
      case _ => false
    }
  }

  /**
   * Returns true if the node is not connected to any edges.
   * @param node Node (set of files) to check.
   * @return true if this set of files is not needed in the graph.
   */
  private def isOrphan(node: QNode) = {
    jobGraph.incomingEdgesOf(node).size == 0 &&
      jobGraph.outgoingEdgesOf(node).size == 0
  }

  /**
   * Utility function for running a method over all function edges.
   * @param edgeFunction Function to run for each FunctionEdge.
   */
  private def foreachFunction(f: (FunctionEdge) => Unit) = {
    jobGraph.edgeSet.toList
      .filter(_.isInstanceOf[FunctionEdge])
      .asInstanceOf[List[FunctionEdge]]
      .sortWith(compare(_,_))
      .foreach(edge => if (running) f(edge))
  }

  private def compare(f1: FunctionEdge, f2: FunctionEdge): Boolean =
    compare(f1.function, f2.function)

  private def compare(f1: QFunction, f2: QFunction): Boolean = {
    val len1 = f1.addOrder.size
    val len2 = f2.addOrder.size
    val len = len1 min len2
    
    for (i <- 0 until len) {
      val order1 = f1.addOrder(i)
      val order2 = f2.addOrder(i)
      if (order1 < order2)
        return true
      if (order1 > order2)
        return false
    }
    if (len1 < len2)
      return true
    else
      return false
  }

  /**
   * Utility function for running a method over all functions, but traversing the nodes in order of dependency.
   * @param edgeFunction Function to run for each FunctionEdge.
   */
  private def traverseFunctions(f: (FunctionEdge) => Unit) {
    val iterator = new TopologicalOrderIterator(this.jobGraph)
    iterator.addTraversalListener(new TraversalListenerAdapter[QNode, QEdge] {
      override def edgeTraversed(event: EdgeTraversalEvent[QNode, QEdge]) = {
        if (running) {
          event.getEdge match {
            case functionEdge: FunctionEdge => f(functionEdge)
            case map: MappingEdge => /* do nothing for mapping functions */
          }
        }
      }
    })
    iterator.foreach(_ => {})
  }

  private def deleteIntermediateOutputs() {
    if (running && !settings.keepIntermediates && success) {
      logger.info("Deleting intermediate files.")
      traverseFunctions(edge => {
        if (edge.function.isIntermediate) {
          logger.debug("Deleting intermediates:" + edge.function.description)
          edge.function.deleteOutputs()
        }
      })
    }
  }

  /**
   * Outputs the graph to a .dot file.
   * http://en.wikipedia.org/wiki/DOT_language
   * @param file Path to output the .dot file.
   */
  private def renderToDot(file: java.io.File) {
    val out = new java.io.FileWriter(file)

    // todo -- we need a nice way to visualize the key pieces of information about commands.  Perhaps a
    // todo -- visualizeString() command, or something that shows inputs / outputs
    val ve = new org.jgrapht.ext.EdgeNameProvider[QEdge] {
      def getEdgeName(function: QEdge) = if (function.dotString == null) "" else function.dotString.replace("\"", "\\\"")
    }

    //val iterator = new TopologicalOrderIterator(qGraph.jobGraph)
    (new DOTExporter(new org.jgrapht.ext.IntegerNameProvider[QNode](), null, ve)).export(out, jobGraph)

    out.close
  }

  /**
   * Returns true if no functions have missing values nor a status of failed.
   * @return true if no functions have missing values nor a status of failed.
   */
  def success = {
    if (numMissingValues > 0) {
      false
    } else if (this.dryRun) {
      true
    } else {
      !this.jobGraph.edgeSet.exists(edge => {
        edge.isInstanceOf[FunctionEdge] && edge.asInstanceOf[FunctionEdge].status == RunnerStatus.FAILED
      })
    }
  }

  def logFailed() {
    foreachFunction(edge => {
      if (edge.status == RunnerStatus.FAILED)
        logEdge(edge)
    })
  }

  /**
   * Returns true if the graph was shutdown instead of exiting on its own.
   */
  def isShutdown = !running

  /**
   * Kills any forked jobs still running.
   */
  def shutdown() {
    // Signal the main thread to shutdown.
    running = false
    // Wait for the thread to finish and exit normally.
    runningLock.synchronized {
      val runners = getRunningJobs.map(_.runner)
      val manager = commandLineManager.asInstanceOf[JobManager[QFunction,JobRunner[QFunction]]]
      if (manager != null) {
        val managerRunners = runners
          .filter(runner => manager.runnerType.isAssignableFrom(runner.getClass))
          .asInstanceOf[List[JobRunner[QFunction]]]
        if (managerRunners.size > 0)
          try {
            manager.tryStop(managerRunners)
          } catch {
            case e => /* ignore */
          }
      }
      runners.foreach(runner =>
        try {
          runner.removeTemporaryFiles()
        } catch {
          case e => /* ignore */
        }
      )
    }
  }
}
