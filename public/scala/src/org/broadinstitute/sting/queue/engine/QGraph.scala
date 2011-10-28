/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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
import org.apache.commons.lang.StringUtils
import org.broadinstitute.sting.queue.util._
import collection.immutable.{TreeSet, TreeMap}
import org.broadinstitute.sting.queue.function.scattergather.{ScatterFunction, CloneFunction, GatherFunction, ScatterGatherableFunction}
import java.util.Date
import org.broadinstitute.sting.utils.Utils
import org.broadinstitute.sting.utils.io.IOUtils

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var settings: QGraphSettings = _

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
  private def managers = List[Any](inProcessManager, commandLineManager)

  private class StatusCounts {
    var pending = 0
    var running = 0
    var failed = 0
    var done = 0
  }
  private val statusCounts = new StatusCounts

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
          val inputs = getQNode(command.inputs.toList.sorted(fileOrdering))
          val outputs = getQNode(command.outputs.toList.sorted(fileOrdering))
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

  private def fillGraph {
    logger.info("Generating graph.")
    fill
    if (settings.dotFile != null)
      renderToDot(settings.dotFile)
    validate()

    if (running && numMissingValues == 0) {
      val scatterGathers = jobGraph.edgeSet.filter(edge => scatterGatherable(edge))
      if (!scatterGathers.isEmpty) {
        logger.info("Generating scatter gather jobs.")

        var addedFunctions = List.empty[QFunction]
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
        fill
        val scatterGatherDotFile = if (settings.expandedDotFile != null) settings.expandedDotFile else settings.dotFile
        if (scatterGatherDotFile != null)
          renderToDot(scatterGatherDotFile)
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
  private def previousFunctions(edge: QEdge): List[FunctionEdge] = {
    var previous = List.empty[FunctionEdge]
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
  private def nextFunctions(edge: QEdge): List[FunctionEdge] = {
    var next = List.empty[FunctionEdge]
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

  private def getReadyJobs(): Set[FunctionEdge] = {
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

    updateGraphStatus(false)

    var readyJobs = getReadyJobs()
    while (running && readyJobs.size > 0) {
      logger.debug("+++++++")
      foreachFunction(readyJobs.toList, edge => {
        if (running) {
          edge.myRunInfo.startTime = new Date()
          edge.getRunInfo.exechosts = Utils.resolveHostname()
          logEdge(edge)
          edge.myRunInfo.doneTime = new Date()
          edge.markAsDone
        }
      })
      readyJobs = getReadyJobs()
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
      if (settings.bsub)
        settings.jobRunner = "Lsf706"
      else if (settings.qsub)
        settings.jobRunner = "GridEngine"
      else if (settings.jobRunner == null)
        settings.jobRunner = "Shell"
      commandLineManager = commandLinePluginManager.createByName(settings.jobRunner)

      for (mgr <- managers) {
        if (mgr != null) {
          val manager = mgr.asInstanceOf[JobManager[QFunction,JobRunner[QFunction]]]
          manager.init()
        }
      }

      if (settings.startFromScratch)
        logger.info("Removing outputs from previous runs.")

      updateGraphStatus(true)

      var readyJobs = TreeSet.empty[FunctionEdge](functionOrdering)
      readyJobs ++= getReadyJobs()
      runningJobs = Set.empty[FunctionEdge]
      var lastRunningCheck = System.currentTimeMillis
      var logNextStatusCounts = true
      var startedJobsToEmail = Set.empty[FunctionEdge]

      while (running && readyJobs.size + runningJobs.size > 0) {

        var startedJobs = Set.empty[FunctionEdge]
        var doneJobs = Set.empty[FunctionEdge]
        var failedJobs = Set.empty[FunctionEdge]

        while (running && readyJobs.size > 0 && !readyRunningCheck(lastRunningCheck)) {
          val edge = readyJobs.head
          edge.runner = newRunner(edge.function)
          edge.start()
          startedJobs += edge
          readyJobs -= edge
          logNextStatusCounts = true
        }

        runningJobs ++= startedJobs
        startedJobsToEmail ++= startedJobs
        statusCounts.pending -= startedJobs.size
        statusCounts.running += startedJobs.size

        if (logNextStatusCounts)
          logStatusCounts
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
          case RunnerStatus.DONE => doneJobs += edge
          case RunnerStatus.FAILED => failedJobs += edge
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

        readyJobs ++= getReadyJobs()
      }

      logStatusCounts
      deleteCleanup(-1)
    } catch {
      case e =>
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

  private def logStatusCounts {
    logger.info("%d Pend, %d Run, %d Fail, %d Done".format(
      statusCounts.pending, statusCounts.running, statusCounts.failed, statusCounts.done))
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
      logger.debug("Already done: " + edge.function.description)
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
   * Checks if the function should have their outptus removed after they finish running
   * @param edges Function to check
   */
  private def addCleanup(edge: FunctionEdge) {
    if (!settings.keepIntermediates)
      if (edge.function.isIntermediate && edge.function.deleteIntermediateOutputs)
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

  private def emailStartedJobs(started: Set[FunctionEdge]) {
    if (settings.statusEmailTo.size > 0) {
      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.subject = "Queue function: Started: " + settings.qSettings.jobNamePrefix
      addStartedFunctions(emailMessage, started.toList)
      emailMessage.trySend(settings.qSettings.emailSettings)
    }
  }

  private def emailFailedJobs(failed: Set[FunctionEdge]) {
    if (settings.statusEmailTo.size > 0) {
      val emailMessage = new EmailMessage
      emailMessage.from = settings.statusEmailFrom
      emailMessage.to = settings.statusEmailTo
      emailMessage.subject = "Queue function: Failure: " + settings.qSettings.jobNamePrefix
      addFailedFunctions(emailMessage, failed.toList)
      emailMessage.trySend(settings.qSettings.emailSettings)
    }
  }

  private def checkRetryJobs(failed: Set[FunctionEdge]) {
    if (settings.retries > 0) {
      for (failedJob <- failed) {
        if (failedJob.function.jobRestartable && failedJob.retries < settings.retries) {
          failedJob.retries += 1
          failedJob.resetToPending(true)
          logger.info("Reset for retry attempt %d of %d: %s".format(
            failedJob.retries, settings.retries, failedJob.function.description))
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

  private def addStartedFunctions(emailMessage: EmailMessage, started: List[FunctionEdge]) {
    if (emailMessage.body == null)
      emailMessage.body = ""
    emailMessage.body += """
    |Started functions:
    |
    |%s
    |""".stripMargin.trim.format(
      started.map(edge => emailDescription(edge)).mkString(nl+nl))
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
      failed.map(edge => emailDescription(edge)).mkString(nl+nl),
      logs.map(_.getAbsolutePath).mkString(nl))

    emailMessage.attachments = logs
  }

  private def emailDescription(edge: FunctionEdge) = {
    val description = new StringBuilder
    if (settings.retries > 0)
      description.append("Attempt %d of %d.%n".format(edge.retries + 1, settings.retries + 1))
    description.append(edge.function.description)
    description.toString()
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

  private def getQNode(files: List[File]) = {
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
  private def foreachFunction(f: (FunctionEdge) => Unit) {
    foreachFunction(jobGraph.edgeSet.toList.filter(_.isInstanceOf[FunctionEdge]).asInstanceOf[List[FunctionEdge]], f)
  }

  /**
   * Utility function for running a method over a list of function edges.
   * @param edegs Edges to traverse.
   * @param edgeFunction Function to run for each FunctionEdge.
   */
  private def foreachFunction(edges: List[FunctionEdge], f: (FunctionEdge) => Unit) {
    edges.sorted(functionOrdering).foreach(edge => if (running) f(edge))
  }

  /**
   * Utility function for running a method over all function edges.
   * @param edgeFunction Function to run for each FunctionEdge.
   */
  private def getFunctionEdges: List[FunctionEdge] = {
    jobGraph.edgeSet.toList.filter(_.isInstanceOf[FunctionEdge]).asInstanceOf[List[FunctionEdge]]
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
            case e => /* ignore */
          }
      }
    }
  }

  /**
   * Returns true if the graph was shutdown instead of exiting on its own.
   */
  def isShutdown = !running

  def getFunctionsAndStatus(functions: List[QFunction]): Map[QFunction, JobRunInfo] = {
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
                case e => /* ignore */
              }
            for (runner <- managerRunners) {
              try {
                runner.cleanup()
              } catch {
                case e => /* ignore */
              }
            }
          } finally {
            try {
              manager.exit()
            } catch {
              case e => /* ignore */
            }
          }
        }
      }
    }
  }
}
