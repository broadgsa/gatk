/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions._
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import org.jgrapht.event.{TraversalListenerAdapter, EdgeTraversalEvent}
import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.queue.function.{InProcessFunction, CommandLineFunction, QFunction}
import org.apache.commons.lang.StringUtils
import org.broadinstitute.gatk.queue.util._
import collection.immutable.{TreeSet, TreeMap}
import org.broadinstitute.gatk.queue.function.scattergather.{ScatterFunction, CloneFunction, GatherFunction, ScatterGatherableFunction}
import java.util.Date
import org.broadinstitute.gatk.utils.Utils
import org.apache.commons.io.{FilenameUtils, FileUtils, IOUtils}
import java.io.{OutputStreamWriter, File}

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var settings: QGraphSettings = _
  var messengers: Seq[QStatusMessenger] = Nil

  private def dryRun = !settings.run
  private var numMissingValues = 0

  private val jobGraph = newGraph
  private val functionOrdering = Ordering.by[FunctionEdge, Iterable[Int]](edge => -graphDepth(edge) +: edge.function.addOrder)
  private val fileOrdering = Ordering.by[File,String](_.getAbsolutePath)
  // A map of nodes by list of files.
  private var nodeMap = TreeMap.empty[Iterable[File], QNode](Ordering.Iterable(fileOrdering))
  // The next unique id for a node if not found in the nodeMap.
  private var nextNodeId = 0

  private var running = true
  private val runningLock = new Object
  private var runningJobs = Set.empty[FunctionEdge]
  private var cleanupJobs = Set.empty[FunctionEdge]

  private val nl = "%n".format()

  private val commandLinePluginManager = new CommandLinePluginManager
  private var commandLineManager: CommandLineJobManager[CommandLineJobRunner] = _
  private val inProcessManager = new InProcessJobManager
  private def managers = Seq[Any](inProcessManager, commandLineManager)

  /**
   * If true, we will write out incremental job reports
   */
  private val INCREMENTAL_JOBS_REPORT = true

  /**
   * Holds the optional jobInfoReporter structure
   */
  private var jobInfoReporter: QJobsReporter = null

  private class StatusCounts {
    var pending = 0
    var running = 0
    var failed = 0
    var done = 0
  }
  private val statusCounts = new StatusCounts

  /**
   * Final initialization step of this QGraph -- tell it runtime setting information
   *
   * The settings aren't necessarily available until after this QGraph object has been constructed, so
   * this function must be called once the QGraphSettings have been filled in.
   *
   * @param settings QGraphSettings
   */
  def initializeWithSettings(settings: QGraphSettings) {
    this.settings = settings
    this.jobInfoReporter = createJobsReporter()
  }

  /**
   * Adds a QScript created CommandLineFunction to the graph.
   * @param command Function to add to the graph.
   */
  def add(command: QFunction) {
    try {
      runningLock.synchronized {
        if (running) {
          command.qSettings = settings.qSettings
          command.freeze()
          val inputFiles = command.inputs
          var outputFiles = command.outputs
          outputFiles :+= command.jobOutputFile
          if (command.jobErrorFile != null)
            outputFiles :+= command.jobErrorFile
          val inputs = getQNode(inputFiles.sorted(fileOrdering))
          val outputs = getQNode(outputFiles.sorted(fileOrdering))
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
        org.broadinstitute.gatk.utils.io.IOUtils.checkTempDir(settings.qSettings.tempDirectory)
        fillGraph()
        val isReady = numMissingValues == 0

        if (this.jobGraph.edgeSet.isEmpty) {
          logger.warn("Nothing to run! Were any Functions added?")
        } else if (settings.getStatus) {
          logger.info("Checking pipeline status.")
          logStatus()
        } else if (this.dryRun) {
          dryRunJobs()
          if (running && isReady) {
            logger.info("Dry run completed successfully!")
            logger.info("Re-run with \"-run\" to execute the functions.")
          }
        } else if (isReady) {
          logger.info("Running jobs.")
          runJobs()
        }

        if (numMissingValues > 0) {
          logger.error("Total missing values: " + numMissingValues)
        }
      }
    }
  }

  private def fillGraph() {
    logger.info("Generating graph.")
    fill()
    if (settings.graphvizFile != null)
      renderGraph(settings.graphvizFile)
    validate()

    if (running && numMissingValues == 0) {
      val scatterGathers = jobGraph.edgeSet.filter(edge => scatterGatherable(edge))
      if (!scatterGathers.isEmpty) {
        logger.info("Generating scatter gather jobs.")

        var addedFunctions = Seq.empty[QFunction]
        for (scatterGather <- scatterGathers) {
          val functions = scatterGather.asInstanceOf[FunctionEdge]
                  .function.asInstanceOf[ScatterGatherableFunction]
                  .generateFunctions()
          addedFunctions ++= functions
        }

        logger.info("Removing original jobs.")
        this.jobGraph.removeAllEdges(scatterGathers)
        prune()

        logger.info("Adding scatter gather jobs.")
        addedFunctions.foreach(function => if (running) this.add(function))

        logger.info("Regenerating graph.")
        fill()
        val scatterGatherDotFile = if (settings.graphvizScatterGatherFile != null) settings.graphvizScatterGatherFile else settings.graphvizFile
        if (scatterGatherDotFile != null)
          renderGraph(scatterGatherDotFile)
        validate()
      }
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
   * Walks up the graph looking for the previous function edges.
   * @param edge Graph edge to examine for the previous functions.
   * @return A list of prior function edges.
   */
  private def previousFunctions(edge: QEdge): Seq[FunctionEdge] = {
    var previous = Seq.empty[FunctionEdge]
    val source = this.jobGraph.getEdgeSource(edge)
    for (incomingEdge <- this.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

        // Stop recursing when we find a function edge and return it
        case functionEdge: FunctionEdge => previous :+= functionEdge

        // For any other type of edge find the jobs preceding the edge
        case edge: QEdge => previous ++= previousFunctions(edge)
      }
    }
    previous
  }

  /**
   * Walks up the graph looking for the next function edges.
   * @param edge Graph edge to examine for the next functions.
   * @return A list of prior function edges.
   */
  private def nextFunctions(edge: QEdge): Seq[FunctionEdge] = {
    var next = Seq.empty[FunctionEdge]
    val target = this.jobGraph.getEdgeTarget(edge)
    for (outgoingEdge <- this.jobGraph.outgoingEdgesOf(target)) {
      outgoingEdge match {

        // Stop recursing when we find a function edge and return it
        case functionEdge: FunctionEdge => next :+= functionEdge

        // For any other type of edge find the jobs following the edge
        case edge: QEdge => next ++= nextFunctions(edge)
      }
    }
    next
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
    asScalaSet(jobGraph.edgeSet).clone().foreach(edge => {
      if (running) edge match {
        case cmd: FunctionEdge => {
          addCollectionOutputs(cmd.outputs)
          addCollectionInputs(cmd.inputs)
        }
        case map: MappingEdge => /* do nothing for mapping edges */
      }
    })
  }

  private def getReadyJobs: Set[FunctionEdge] = {
    jobGraph.edgeSet.filter{
      case f: FunctionEdge =>
        this.previousFunctions(f).forall(_.status == RunnerStatus.DONE) && f.status == RunnerStatus.PENDING
      case _ => false
    }.toSet.asInstanceOf[Set[FunctionEdge]]
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

    if (running) {
      for (orphan <- jobGraph.vertexSet.filter(isOrphan(_))) {
        jobGraph.removeVertex(orphan)
        nodeMap -= orphan.files
      }
    }
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
    if (settings.startFromScratch)
      logger.info("Will remove outputs from previous runs.")

    updateGraphStatus(cleanOutputs = false)

    var readyJobs = getReadyJobs
    while (running && readyJobs.size > 0) {
      logger.debug("+++++++")
      foreachFunction(readyJobs.toSeq, edge => {
        if (running) {
          edge.myRunInfo.startTime = new Date()
          edge.getRunInfo.exechosts = Utils.resolveHostname()
          logEdge(edge)
          edge.myRunInfo.doneTime = new Date()
          edge.markAsDone()
        }
      })
      readyJobs = getReadyJobs
    }
  }

  private def logEdge(edge: FunctionEdge) {
    logger.info("-------")
    logger.info("%-8s %s".format(StringUtils.capitalize(edge.status.toString) + ":", edge.function.description))
    if (logger.isDebugEnabled) {
      logger.debug("Inputs:  " + edge.inputs)
      logger.debug("Outputs: " + edge.outputs)
      logger.debug("Done+:   " + edge.function.doneOutputs.filter(_.exists()))
      logger.debug("Done-:   " + edge.function.doneOutputs.filterNot(_.exists()))
      logger.debug("CmdDir:  " + edge.function.commandDirectory)
      logger.debug("Temp?:   " + edge.function.isIntermediate)
      logger.debug("Prev:    " +
        (if (edge.resetFromStatus == null) "none" else StringUtils.capitalize(edge.resetFromStatus.toString)) +
        " (reset = " + (edge.resetFromStatus != null && edge.resetFromStatus != edge.status) + ")" )
    }
    logger.info("Log:     " + edge.function.jobOutputFile.getAbsolutePath)
    if (edge.function.jobErrorFile != null)
      logger.info("Error:   " + edge.function.jobErrorFile.getAbsolutePath)
  }

  /**
   * Logs job statuses by traversing the graph and looking for status-related files
   */
  private def logStatus() {
    updateGraphStatus(cleanOutputs = false)
    doStatus(status => logger.info(status))
  }

  /**
   * Runs the jobs by traversing the graph.
   */
  private def runJobs() {
    try {
      if (settings.bsub) {
        settings.jobRunner = "Lsf706"
      }
      else if (settings.qsub || settings.qsubBroad) {
        settings.jobRunner = "GridEngine"
        if ( settings.qsubBroad ) {
          settings.qSettings.useBroadClusterSettings = true
        }
      }
      else if (settings.jobRunner == null) {
        settings.jobRunner = "Shell"
      }
      commandLineManager = commandLinePluginManager.createByName(settings.jobRunner)

      for (mgr <- managers) {
        if (mgr != null) {
          val manager = mgr.asInstanceOf[JobManager[QFunction,JobRunner[QFunction]]]
          manager.init()
        }
      }

      if (settings.startFromScratch)
        logger.info("Removing outputs from previous runs.")

      updateGraphStatus(cleanOutputs = true)

      var readyJobs = TreeSet.empty[FunctionEdge](functionOrdering)
      readyJobs ++= getReadyJobs
      runningJobs = Set.empty[FunctionEdge]
      var lastRunningCheck = System.currentTimeMillis
      var logNextStatusCounts = true
      var startedJobsToEmail = Set.empty[FunctionEdge]

      while (running && readyJobs.size + runningJobs.size > 0) {

        var startedJobs = Set.empty[FunctionEdge]
        var doneJobs = Set.empty[FunctionEdge]
        var failedJobs = Set.empty[FunctionEdge]

        def startJobs: Boolean = {

          def canRunMoreConcurrentJobs: Boolean =
            if(settings.maximumNumberOfConcurrentJobs.isDefined)
              runningJobs.size + startedJobs.size < settings.maximumNumberOfConcurrentJobs.get
            else
              true

          running && readyJobs.size > 0 &&
            !readyRunningCheck(lastRunningCheck) &&
            canRunMoreConcurrentJobs
        }

        while (startJobs) {
          val edge = readyJobs.head
          edge.runner = newRunner(edge.function)
          edge.start()
          messengers.foreach(_.started(jobShortName(edge.function)))
          startedJobs += edge
          readyJobs -= edge
          logNextStatusCounts = true
        }

        runningJobs ++= startedJobs
        startedJobsToEmail ++= startedJobs
        statusCounts.pending -= startedJobs.size
        statusCounts.running += startedJobs.size

        if (logNextStatusCounts)
          logStatusCounts()
        logNextStatusCounts = false

        deleteCleanup(lastRunningCheck)

        if (running && startedJobs.size > 0 && !readyRunningCheck(lastRunningCheck)) {
          emailStartedJobs(startedJobsToEmail)
          startedJobsToEmail = Set.empty[FunctionEdge]
        }

        if (readyJobs.size == 0 && runningJobs.size > 0) {
          runningLock.synchronized {
            if (running) {
              val timeout = nextRunningCheck(lastRunningCheck)
              if (timeout > 0)
                runningLock.wait(timeout)
            }
          }
        }

        lastRunningCheck = System.currentTimeMillis
        updateStatus()

        runningJobs.foreach(edge => edge.status match {
          case RunnerStatus.DONE => {
            doneJobs += edge
            messengers.foreach(_.done(jobShortName(edge.function)))
          }
          case RunnerStatus.FAILED => {
            failedJobs += edge
            messengers.foreach(_.exit(jobShortName(edge.function), edge.function.jobErrorLines.mkString("%n".format())))
          }
          case RunnerStatus.RUNNING => /* do nothing while still running */
        })

        runningJobs --= doneJobs
        runningJobs --= failedJobs

        startedJobsToEmail &~= failedJobs

        addCleanup(doneJobs)

        statusCounts.running -= doneJobs.size
        statusCounts.running -= failedJobs.size
        statusCounts.done += doneJobs.size
        statusCounts.failed += failedJobs.size

        if (doneJobs.size > 0 || failedJobs.size > 0)
          logNextStatusCounts = true

        if (running && failedJobs.size > 0) {
          emailFailedJobs(failedJobs)
          checkRetryJobs(failedJobs)
        }

        // incremental
        if ( logNextStatusCounts && INCREMENTAL_JOBS_REPORT ) {
          logger.info("Writing incremental jobs reports...")
          writeJobsReport(plot = false)
        }

        readyJobs ++= getReadyJobs
      }

      logStatusCounts()
      deleteCleanup(-1)
    } catch {
      case e: Throwable =>
        logger.error("Uncaught error running jobs.", e)
        throw e
    } finally {
      emailStatus()
    }
  }

  private def readyRunningCheck(lastRunningCheck: Long) =
    lastRunningCheck > 0 && nextRunningCheck(lastRunningCheck) <= 0

  private def nextRunningCheck(lastRunningCheck: Long) =
    ((30 * 1000L) - (System.currentTimeMillis - lastRunningCheck))

  def formattedStatusCounts: String = {
    "%d Pend, %d Run, %d Fail, %d Done".format(
      statusCounts.pending, statusCounts.running, statusCounts.failed, statusCounts.done)
  }

  private def logStatusCounts() {
    logger.info(formattedStatusCounts)
  }

  /**
   * Updates the status of edges in the graph.
   * @param cleanOutputs If true will delete outputs when setting edges to pending.
   */
  private def updateGraphStatus(cleanOutputs: Boolean) {
    if (settings.startFromScratch)
      foreachFunction(edge => edge.resetToPending(cleanOutputs))
    else
      traverseFunctions(edge => checkDone(edge, cleanOutputs))
    traverseFunctions(edge => recheckDone(edge))
  }

  // TODO: Yet another field to add (with overloads) to QFunction?
  private def jobShortName(function: QFunction): String = {
    var name = function.analysisName
    if (function.isInstanceOf[CloneFunction]) {
      val cloneFunction = function.asInstanceOf[CloneFunction]
      name += " %d of %d".format(cloneFunction.cloneIndex, cloneFunction.cloneCount)
    }
    name
  }

  /**
   * First pass that checks if an edge is done or if it's an intermediate edge if it can be skipped.
   * This function may modify the status of previous edges if it discovers that the edge passed in
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
   * Second pass which
   * a) Updates the status counts based on the function statuses
   * b) Checks if the edge is a completed intermediate edge then adds it to the set of candidates for cleanup
   * @param edge Edge to check to see if it's done or skipped.
   */
  private def recheckDone(edge: FunctionEdge) {
    edge.status match {
      case RunnerStatus.PENDING => statusCounts.pending += 1
      case RunnerStatus.FAILED => statusCounts.failed += 1
      case RunnerStatus.DONE => statusCounts.done += 1
      case RunnerStatus.SKIPPED => statusCounts.done += 1
    }
    
    if (edge.status == RunnerStatus.DONE || edge.status == RunnerStatus.SKIPPED) {
      if (logger.isDebugEnabled)
        logEdge(edge)
      addCleanup(edge)
    }
  }

  /**
   * Checks if the functions should have their outptus removed after they finish running
   * @param edges Functions to check
   */
  private def addCleanup(edges: Traversable[FunctionEdge]) {
    edges.foreach(addCleanup(_))
  }

  /**
   * Checks if the function should have their outputs removed after they finish running
   * @param edge Function to check
   */
  private def addCleanup(edge: FunctionEdge) {
    if (!settings.keepIntermediates)
      if (edge.function.isIntermediate)
        cleanupJobs += edge
  }

  /**
   * Continues deleting the outputs of intermediate jobs that are no longer needed until it's time to recheck running status.
   * @param lastRunningCheck The last time the status was checked.
   */
  private def deleteCleanup(lastRunningCheck: Long) {
    var doneJobs = Set.empty[FunctionEdge]

    for (edge <- cleanupJobs) {
      val nextDone = nextFunctions(edge).forall(next => {
        val status = next.status
        (status == RunnerStatus.DONE || status == RunnerStatus.SKIPPED)
      })

      if (nextDone)
        doneJobs += edge
    }

    for (edge <- doneJobs) {
      if (running && !readyRunningCheck(lastRunningCheck)) {
        logger.debug("Deleting intermediates:" + edge.function.description)
        edge.function.deleteOutputs()
        cleanupJobs -= edge
      }
    }
  }

  /**
   * Returns the graph depth for the function.
   * @param edge Function edge to get the edge for.
   * @return the graph depth for the function.
   */
  private def graphDepth(edge: FunctionEdge): Int = {
    if (edge.depth < 0) {
      val previous = previousFunctions(edge)
      if (previous.size == 0)
        edge.depth = 0
      else
        edge.depth = previous.map(f => graphDepth(f)).max + 1
    }
    edge.depth
  }

  /**
   * From the previous edges, resets any that are marked as skipped to pending.
   * If those that are reset have skipped edges, those skipped edges are recursively also set
   * to pending.
   * Any edges after this edge are also reset to pending.
   * @param edge Dependent edge.
   * @param previous Previous edges that provide inputs to edge.
   * @param cleanOutputs If true will clean up the output files when resetting jobs to pending.
   */
  private def resetPreviousSkipped(edge: FunctionEdge, previous: Seq[FunctionEdge], cleanOutputs: Boolean) {
    val edges = previous.filter(_.status == RunnerStatus.SKIPPED) ++ this.nextFunctions(edge).filter(_.status != RunnerStatus.PENDING)
    for (resetEdge <- edges) {
      resetEdge.resetToPending(cleanOutputs)
      resetPreviousSkipped(resetEdge, this.previousFunctions(resetEdge), cleanOutputs)
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

  private def emailStartedJobs(started: Set[FunctionEdge]) {
    if (settings.statusEmailTo.size > 0) {
      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.subject = "Queue function: Started: " + settings.qSettings.runName
      addStartedFunctions(emailMessage, started.toSeq)
      emailMessage.trySend(settings.emailSettings)
    }
  }

  private def emailFailedJobs(failed: Set[FunctionEdge]) {
    if (settings.statusEmailTo.size > 0) {
      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.subject = "Queue function: Failure: " + settings.qSettings.runName
      addFailedFunctions(emailMessage, failed.toSeq)
      emailMessage.trySend(settings.emailSettings)
    }
  }

  private def checkRetryJobs(failed: Set[FunctionEdge]) {
    if (settings.retries > 0) {
      for (failedJob <- failed) {
        if (failedJob.function.jobRestartable && failedJob.function.retries < settings.retries) {
          failedJob.function.retries += 1
          failedJob.function.setupRetry()
          failedJob.resetToPending(cleanOutputs = true)
          logger.info("Reset for retry attempt %d of %d: %s".format(
            failedJob.function.retries, settings.retries, failedJob.function.description))
          statusCounts.failed -= 1
          statusCounts.pending += 1
        } else {
          logger.info("Giving up after retrying %d times: %s".format(
            settings.retries, failedJob.function.description))
        }
      }
    }
  }

  private def emailStatus() {
    if (running && settings.statusEmailTo.size > 0) {
      var failed = Seq.empty[FunctionEdge]
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
        emailMessage.subject = "Queue run: Success: " + settings.qSettings.runName
      } else {
        emailMessage.subject = "Queue run: Failure: " + settings.qSettings.runName
        addFailedFunctions(emailMessage, failed)
      }
      emailMessage.trySend(settings.emailSettings)
    }
  }

  private def addStartedFunctions(emailMessage: EmailMessage, started: Seq[FunctionEdge]) {
    if (emailMessage.body == null)
      emailMessage.body = ""
    emailMessage.body += """
    |Started functions:
    |
    |%s
    |""".stripMargin.trim.format(
      started.map(edge => emailDescription(edge)).mkString(nl+nl))
  }

  private def addFailedFunctions(emailMessage: EmailMessage, failed: Seq[FunctionEdge]) {
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
      failed.map(edge => emailDescription(edge)).mkString(nl+nl),
      logs.map(_.getAbsolutePath).mkString(nl))

    emailMessage.attachments = logs
  }

  private def emailDescription(edge: FunctionEdge) = {
    val description = new StringBuilder
    if (settings.retries > 0)
      description.append("Attempt %d of %d.%n".format(edge.function.retries + 1, settings.retries + 1))
    description.append(edge.function.description)
    description.toString()
  }

  private def logFiles(edge: FunctionEdge) = {
    var failedOutputs = Seq.empty[File]
    failedOutputs :+= edge.function.jobOutputFile
    if (edge.function.jobErrorFile != null)
      failedOutputs :+= edge.function.jobErrorFile
    failedOutputs.filter(file => file != null && file.exists)
  }

  /**
   * Tracks analysis status.
   */
  private class AnalysisStatus(val analysisName: String) {
    val jobs = new GroupStatus
    val scatter = new GroupStatus
    val gather = new GroupStatus

    def total = jobs.total + scatter.total + gather.total
    def done = jobs.done + scatter.done + gather.done
    def failed = jobs.failed + scatter.failed + gather.failed
    def skipped = jobs.skipped + scatter.skipped + gather.skipped
  }

  /**
   * Tracks status of a group of jobs.
   */
  private class GroupStatus {
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
    buffer.toString()
  }

  /**
   * Gets job statuses by traversing the graph and looking for status-related files
   */
  private def doStatus(statusFunc: String => Unit) {
    var statuses = Seq.empty[AnalysisStatus]
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
      val total = status.total
      val done = status.done
      val failed = status.failed
      val skipped = status.skipped
      val jobsTotal = status.jobs.total
      val jobsDone = status.jobs.done
      val gatherTotal = status.gather.total
      val gatherDone = status.gather.done

      var summaryStatus = RunnerStatus.PENDING
      if (failed > 0)
        summaryStatus = RunnerStatus.FAILED
      else if (gatherDone == gatherTotal && jobsDone == jobsTotal)
        summaryStatus = RunnerStatus.DONE
      else if (done + skipped == total)
        summaryStatus = RunnerStatus.SKIPPED
      else if (done > 0)
        summaryStatus = RunnerStatus.RUNNING

      var info = ("%-" + maxWidth + "s %7s")
              .format(status.analysisName, "[" + summaryStatus.toString + "]")
      if (status.jobs.total > 1) {
        info += formatGroupStatus(status.jobs)
      }
      if (status.scatter.total + status.gather.total > 1) {
        info += formatGroupStatus(status.scatter, "s:")
        info += formatGroupStatus(status.gather, "g:")
      }
      statusFunc(info)
    })
  }

  /**
   * Updates a status map with scatter/gather status information (e.g. counts)
   */
  private def updateAnalysisStatus(stats: AnalysisStatus, edge: FunctionEdge) {
    if (edge.function.isInstanceOf[ScatterFunction]) {
      updateGroupStatus(stats.scatter, edge)
    } else if (edge.function.isInstanceOf[CloneFunction]) {
      updateGroupStatus(stats.scatter, edge)
    } else if (edge.function.isInstanceOf[GatherFunction]) {
      updateGroupStatus(stats.gather, edge)
    } else {
      updateGroupStatus(stats.jobs, edge)
    }
  }

  private def updateGroupStatus(groupStatus: GroupStatus, edge: FunctionEdge) {
    groupStatus.total += 1
    edge.status match {
      case RunnerStatus.DONE => groupStatus.done += 1
      case RunnerStatus.FAILED => groupStatus.failed += 1
      case RunnerStatus.SKIPPED => groupStatus.skipped += 1
      /* can't tell the difference between pending and running right now! */
      case RunnerStatus.PENDING =>
      case RunnerStatus.RUNNING =>
    }
  }

  /**
   * Formats a status into nice strings
   */
  private def formatGroupStatus(stats: GroupStatus, prefix: String = "") = {
    " %s%dt/%dd/%df".format(
      prefix, stats.total, stats.done, stats.failed)
  }

  /**
   *   Creates a new graph where if new edges are needed (for cyclic dependency checking) they can be automatically created using a generic MappingFunction.
   * @return A new graph
   */
  private def newGraph = new SimpleDirectedGraph[QNode, QEdge](new EdgeFactory[QNode, QEdge] {
    def createEdge(input: QNode, output: QNode) = new MappingEdge(input, output)})

  private def getQNode(files: Seq[File]) = {
    nodeMap.get(files) match {
      case Some(node) =>
        node
      case None =>
        if (nextNodeId % 100 == 0)
          logger.debug("adding QNode: " + nextNodeId)
        val node = new QNode(nextNodeId, files)
        nextNodeId += 1
        jobGraph.addVertex(node)
        nodeMap += files -> node
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
          val input = getQNode(Seq(file))
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
          val output = getQNode(Seq(file))
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
   * @param f Function to run for each FunctionEdge.
   */
  private def foreachFunction(f: (FunctionEdge) => Unit) {
    foreachFunction(jobGraph.edgeSet.toSeq.filter(_.isInstanceOf[FunctionEdge]).asInstanceOf[Seq[FunctionEdge]], f)
  }

  /**
   * Utility function for running a method over a list of function edges.
   * @param edges Edges to traverse.
   * @param f Function to run for each FunctionEdge.
   */
  private def foreachFunction(edges: Seq[FunctionEdge], f: (FunctionEdge) => Unit) {
    edges.sorted(functionOrdering).foreach(edge => if (running) f(edge))
  }

  /**
   * Utility function returning all function edges.
   */
  private def getFunctionEdges: Seq[FunctionEdge] = {
    jobGraph.edgeSet.toSeq.filter(_.isInstanceOf[FunctionEdge]).asInstanceOf[Seq[FunctionEdge]]
  }

  /**
   * Utility function for running a method over all functions, but traversing the nodes in order of dependency.
   * @param f Function to run for each FunctionEdge.
   */
  private def traverseFunctions(f: (FunctionEdge) => Unit) {
    val iterator = new TopologicalOrderIterator(this.jobGraph)
    iterator.addTraversalListener(new TraversalListenerAdapter[QNode, QEdge] {
      override def edgeTraversed(event: EdgeTraversalEvent[QNode, QEdge]) {
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

  /**
   * Outputs the graph to a .gv DOT file.
   * http://www.graphviz.org/Documentation.php
   * http://en.wikipedia.org/wiki/DOT_language
   * @param file Path to output the .gv file.
   */
  private def renderGraph(file: java.io.File) {
    val vertexIDProvider = new org.jgrapht.ext.VertexNameProvider[QNode] {
      def getVertexName(node: QNode) = node.id.toString
    }

    val vertexLabelProvider = new org.jgrapht.ext.VertexNameProvider[QNode] {
      // The QGraph fills in with single file nodes between nodes that contain more than one file.
      // We only need to display the single element nodes.
      def getVertexName(node: QNode) = {
        if (!node.files.isEmpty && node.files.tail.isEmpty)
          node.files.head.getName
        else
          ""
      }
    }

    val edgeNameProvider = new org.jgrapht.ext.EdgeNameProvider[QEdge] {
      def getEdgeName(edge: QEdge) = {
        if (edge.shortDescription != null)
          edge.shortDescription.replace("\"", "\\\"")
        else
          ""
      }
    }

    val exporter = new DOTExporter(vertexIDProvider, vertexLabelProvider, edgeNameProvider)

    val out = new OutputStreamWriter(FileUtils.openOutputStream(file))
    try {
      exporter.export(out, jobGraph)
    } finally {
      IOUtils.closeQuietly(out)
    }
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
        if (edge.isInstanceOf[FunctionEdge]) {
          val status = edge.asInstanceOf[FunctionEdge].status
          (status == RunnerStatus.PENDING || status == RunnerStatus.RUNNING || status == RunnerStatus.FAILED)
        } else {
          false
        }
      })
    }
  }

  def logFailed() {
    foreachFunction(edge => {
      if (edge.status == RunnerStatus.FAILED)
        logEdge(edge)
    })
  }


  private def updateStatus() {
    val runners = runningJobs.map(_.runner)
    for (mgr <- managers) {
      if (mgr != null) {
        val manager = mgr.asInstanceOf[JobManager[QFunction,JobRunner[QFunction]]]
        val managerRunners = runners
          .filter(runner => manager.runnerType.isAssignableFrom(runner.getClass))
          .asInstanceOf[Set[JobRunner[QFunction]]]
        if (managerRunners.size > 0)
          try {
            val updatedRunners = manager.updateStatus(managerRunners)
            for (runner <- managerRunners.diff(updatedRunners)) {
              runner.checkUnknownStatus()
            }
          } catch {
            case e: Throwable => /* ignore */
          }
      }
    }
  }

  /**
   * Create the jobsReporter for this QGraph, based on the settings data.
   *
   * Must be called after settings has been initialized properly
   *
   * @return
   */
  private def createJobsReporter(): QJobsReporter = {
    val jobStringName = if (settings.jobReportFile != null)
      settings.jobReportFile
    else
      settings.qSettings.runName + ".jobreport.txt"

    val reportFile = org.broadinstitute.gatk.utils.io.IOUtils.absolute(settings.qSettings.runDirectory, jobStringName)

    val pdfFile = if ( settings.run )
      Some(org.broadinstitute.gatk.utils.io.IOUtils.absolute(settings.qSettings.runDirectory, FilenameUtils.removeExtension(jobStringName) + ".pdf"))
    else
      None

    new QJobsReporter(settings.disableJobReport, reportFile, pdfFile)
  }

  /**
   * Write, if possible, the jobs report
   */
  def writeJobsReport(plot: Boolean = true) {
    // note: the previous logic didn't write the job report if the system was shutting down, but I don't
    // see any reason not to write the job report
    if ( jobInfoReporter != null )
      jobInfoReporter.write(this, plot)
  }

  /**
   * Returns true if the graph was shutdown instead of exiting on its own.
   */
  def isShutdown = !running

  def getFunctionsAndStatus: Map[QFunction, JobRunInfo] = {
    getFunctionEdges.map(edge => (edge.function, edge.getRunInfo)).toMap
  }

  /**
   * Kills any forked jobs still running.
   */
  def shutdown() {
    // Signal the main thread to shutdown.
    running = false

    // Try and wait for the thread to finish and exit normally.
    runningLock.synchronized {
      runningLock.notify()
    }

    // Start killing jobs.
    runningLock.synchronized {
      val runners = runningJobs.map(_.runner)
      runningJobs = Set.empty[FunctionEdge]
      for (mgr <- managers) {
        if (mgr != null) {
          val manager = mgr.asInstanceOf[JobManager[QFunction,JobRunner[QFunction]]]
          try {
            val managerRunners = runners
              .filter(runner => manager.runnerType.isAssignableFrom(runner.getClass))
              .asInstanceOf[Set[JobRunner[QFunction]]]
            if (managerRunners.size > 0)
              try {
                manager.tryStop(managerRunners)
              } catch {
                case e: Throwable => /* ignore */
              }
            for (runner <- managerRunners) {
              try {
                runner.cleanup()
              } catch {
                case e: Throwable => /* ignore */
              }
            }
          } finally {
            try {
              manager.exit()
            } catch {
              case e: Throwable => /* ignore */
            }
          }
        }
      }
    }
  }
}
