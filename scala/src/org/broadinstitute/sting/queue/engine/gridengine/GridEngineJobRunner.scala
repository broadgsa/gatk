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

package org.broadinstitute.sting.queue.engine.gridengine

import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.queue.util.{Logging,Retry}
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.engine.{RunnerStatus, CommandLineJobRunner}
import org.ggf.drmaa.{DrmaaException,JobInfo,JobTemplate,Session,SessionFactory}
import java.util.Collections

/**
 * Runs jobs on a Grid Engine compute cluster.
 */
class GridEngineJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {
  // Run the static initializer for GridEngineJobRunner
  GridEngineJobRunner

  /** Job Id of the currently executing job. */
  private var jobId: String = _

  /** Last known status */
  private var lastStatus: RunnerStatus.Value = _

  /** The last time the status was updated */
  protected var lastStatusUpdate: Long = _

  def start() {
    GridEngineJobRunner.gridEngineSession.synchronized {
      val gridEngineJob: JobTemplate = GridEngineJobRunner.gridEngineSession.createJobTemplate

      // Force the remote environment to inherit local environment settings
      var nativeSpecString: String = "-V"

      // Set the display name to < 512 characters of the description
      // NOTE: Not sure if this is configuration specific?
      gridEngineJob.setJobName(GridEngineJobRunner.toJobName(function.description.take(500)))

      // Set the output file for stdout
      gridEngineJob.setOutputPath(":" + function.jobOutputFile.getPath)

      // Set the current working directory
      gridEngineJob.setWorkingDirectory(function.commandDirectory.getPath)

      // If the error file is set specify the separate output for stderr
      // Otherwise join with stdout
      if (Option(function.jobErrorFile) != None) {
        gridEngineJob.setErrorPath(":" + function.jobErrorFile.getPath)
      } else {
        gridEngineJob.setJoinFiles(true)
      }

      // If a project name is set specify the project name
      if (Option(function.jobProject) != None) {
        nativeSpecString += " -P " + function.jobProject
      }

      // If the job queue is set specify the job queue
      if (Option(function.jobQueue) != None) {
        nativeSpecString += " -q " + function.jobQueue
      }

      // If the memory limit is set (GB) specify the memory limit
      if (function.memoryLimit.isDefined) {
        val memAvl: String = function.memoryLimit.get + "G"
        val memMax: String = (function.memoryLimit.get * 1.2 * 1024).ceil.toInt + "M"
        nativeSpecString += " -l mem_free=" + memAvl + ",h_rss=" + memMax
      }

      // If the priority is set (user specified Int) specify the priority
      if (function.jobPriority.isDefined) {
        nativeSpecString += " -p " + function.jobPriority.get
      }

      gridEngineJob.setNativeSpecification(nativeSpecString)

      // Instead of running the function.commandLine, run "sh <jobScript>"
      gridEngineJob.setRemoteCommand("sh")
      gridEngineJob.setArgs(Collections.singletonList(jobScript.toString))

      // Allow advanced users to update the request via QFunction.updateJobRun()
      updateJobRun(gridEngineJob)

      updateStatus(RunnerStatus.RUNNING)

      // Start the job and store the id so it can be killed in tryStop
      try {
        Retry.attempt(() => {
          try {
            jobId = GridEngineJobRunner.gridEngineSession.runJob(gridEngineJob)
          } catch {
            case de: DrmaaException => throw new QException("Unable to submit job: " + de.getLocalizedMessage)
          }
        }, 1, 5, 10)
      } finally {
        // Prevent memory leaks
        GridEngineJobRunner.gridEngineSession.deleteJobTemplate(gridEngineJob)
      }
      logger.info("Submitted Grid Engine job id: " + jobId)
    }
  }

  def status = this.lastStatus

  private def updateStatus(updatedStatus: RunnerStatus.Value) {
    this.lastStatus = updatedStatus
    this.lastStatusUpdate = System.currentTimeMillis
  }
}

object GridEngineJobRunner extends Logging {
  private val gridEngineSession = SessionFactory.getFactory.getSession

  /** Amount of time a job can go without status before giving up. */
  private val unknownStatusMaxSeconds = 5 * 60

  initGridEngine()

  /**
   * Initialize the Grid Engine library.
   */
  private def initGridEngine() {
    gridEngineSession.synchronized {
      try {
        gridEngineSession.init("")
      } catch {
        case de: DrmaaException =>
          logger.error("Issue initializing Grid Engine", de)
          throw new QException("init() failed", de)
      }
    }
  }

  /**
   * Updates the status of a list of jobs.
   * @param runners Runners to update.
   */
  def updateStatus(runners: Set[GridEngineJobRunner]) {
    var updatedRunners = Set.empty[GridEngineJobRunner]
    gridEngineSession.synchronized {
      runners.foreach(runner => if (updateRunnerStatus(runner)) {updatedRunners += runner})
    }

    for (runner <- runners.diff(updatedRunners)) {
      checkUnknownStatus(runner)
    }
  }

  /**
   * Tries to stop any running jobs.
   * @param runners Runners to stop.
   */
  def tryStop(runners: Set[GridEngineJobRunner]) {
    // Stop runners. SIGTERM(15) is preferred to SIGKILL(9).
    // Only way to send SIGTERM is for the Sys Admin set the terminate_method
    // resource of the designated queue to SIGTERM
    gridEngineSession.synchronized {
      for (runner <- runners.filterNot(runner => Option(runner.jobId) == None)) {
        try {
          gridEngineSession.control(runner.jobId, Session.TERMINATE)
        } catch {
          case e =>
            logger.error("Unable to kill job " + runner.jobId, e)
        }
      }
      gridEngineSession.exit()
    }
  }

  private def updateRunnerStatus(runner: GridEngineJobRunner): Boolean = {
    var returnStatus: RunnerStatus.Value = null

    try {
      val jobStatus = gridEngineSession.getJobProgramStatus(runner.jobId);
      jobStatus match {
        case Session.QUEUED_ACTIVE => returnStatus = RunnerStatus.RUNNING
        case Session.DONE =>
          val jobInfo: JobInfo = gridEngineSession.wait(runner.jobId, Session.TIMEOUT_NO_WAIT)
          if ((jobInfo.hasExited && jobInfo.getExitStatus > 0)
              || jobInfo.hasSignaled
              || jobInfo.wasAborted)
            returnStatus = RunnerStatus.FAILED
          else
            returnStatus = RunnerStatus.DONE
        case Session.FAILED => returnStatus = RunnerStatus.FAILED
        case Session.UNDETERMINED => logger.warn("Unable to determine status of Grid Engine job id " + runner.jobId)
        case _ => returnStatus = RunnerStatus.RUNNING
      }
    } catch {
      // getJobProgramStatus will throw an exception once wait has run, as the
      // job will be reaped.  If the status is currently DONE or FAILED, return
      // the status.
      case de: DrmaaException =>
        if (runner.lastStatus == RunnerStatus.DONE || runner.lastStatus == RunnerStatus.FAILED)
          returnStatus = runner.lastStatus
        else
          logger.warn("Unable to determine status of Grid Engine job id " + runner.jobId, de)
    }

    Option(returnStatus) match {
      case Some(returnStatus) =>
        runner.updateStatus(returnStatus)
        return true
      case None => return false
    }
  }

  private def checkUnknownStatus(runner: GridEngineJobRunner) {
    val unknownStatusSeconds = (System.currentTimeMillis - runner.lastStatusUpdate)
    if (unknownStatusSeconds > (unknownStatusMaxSeconds * 1000L)) {
      // Unknown status has been returned for a while now.
      runner.updateStatus(RunnerStatus.FAILED)
      logger.error("Unable to read Grid Engine status for %d minutes: job id %d: %s".format(unknownStatusSeconds/60, runner.jobId, runner.function.description))
    }
  }

  // Reap what we've sown
  override def finalize() {
    gridEngineSession.exit()
  }

  // Grid Engine disallows certain characters from being in job names.
  // This replaces all illegal characters with underscores
  private def toJobName(name: String): String = {
    name.replaceAll("""[\n\t\r/:@\\*?]""", "_")
  }
}
