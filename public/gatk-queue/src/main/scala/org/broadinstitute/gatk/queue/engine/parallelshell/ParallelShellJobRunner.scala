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

package org.broadinstitute.gatk.queue.engine.parallelshell

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.{ RunnerStatus, CommandLineJobRunner }
import java.util.Date
import org.broadinstitute.gatk.utils.Utils
import org.broadinstitute.gatk.utils.runtime.{ ProcessSettings, OutputStreamSettings }
import scala.concurrent._
import ExecutionContext.Implicits.global
import scala.util.{ Success, Failure }
import org.broadinstitute.gatk.queue.util.Logging

/**
 * Runs multiple jobs locally without blocking.
 * Use this with care as it might not be the most efficient way to run things.
 * However, for some scenarios, such as running multiple single threaded
 * programs concurrently it can be quite useful.
 *
 * All this code is based on the normal shell runner in GATK Queue and all
 * credits for everything except the concurrency part goes to the GATK team.
 *
 * @author Johan Dahlberg - 20150611
 *
 * @param function Command to run.
 */
class ParallelShellJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {

  // Controller on the thread that started the job
  val controller: ThreadSafeProcessController = new ThreadSafeProcessController()

  // Once the application exits this promise will be fulfilled.
  val finalExitStatus = Promise[Int]()

  /**
   * Runs the function on the local shell.
   */
  def start() {
    val commandLine = Array("sh", jobScript.getAbsolutePath)
    val stdoutSettings = new OutputStreamSettings
    val stderrSettings = new OutputStreamSettings
    val mergeError = function.jobErrorFile == null

    stdoutSettings.setOutputFile(function.jobOutputFile, true)
    if (function.jobErrorFile != null)
      stderrSettings.setOutputFile(function.jobErrorFile, true)

    if (logger.isDebugEnabled) {
      stdoutSettings.printStandard(true)
      stderrSettings.printStandard(true)
    }

    val processSettings = new ProcessSettings(
      commandLine, mergeError, function.commandDirectory, null,
      null, stdoutSettings, stderrSettings)

    updateJobRun(processSettings)

    getRunInfo.startTime = new Date()
    getRunInfo.exechosts = Utils.resolveHostname()
    updateStatus(RunnerStatus.RUNNING)

    // Run the command line process in a future.
    val executedFuture =
      future { controller.exec(processSettings) }

    // Register a callback on the completion of the future, making sure that
    // the status of the job is updated accordingly. 
    executedFuture.onComplete {
      case Success(exitStatus) =>
        logger.debug(commandLine.mkString(" ") + " :: Got return on exit status in future: " + exitStatus)
        finalExitStatus.success(exitStatus)
        getRunInfo.doneTime = new Date()
        exitStatusUpdateJobRunnerStatus(exitStatus)
      case Failure(throwable) =>
        logger.debug(
          "Failed in return from run with: " +
            throwable.getClass.getCanonicalName + " :: " +
            throwable.getMessage)
        finalExitStatus.failure(throwable)
        getRunInfo.doneTime = new Date()
        updateStatus(RunnerStatus.FAILED)
    }
  }

  /**
   * Possibly invoked from a shutdown thread, find and
   * stop the controller from the originating thread
   */
  def tryStop() = {
    try {
      controller.tryDestroy()
    } catch {
      case e: Exception =>
        logger.error("Unable to kill shell job: " + function.description, e)
    }
  }

  /**
   * Update the status of the runner based on the exit status
   * of the process.
   */
  def exitStatusUpdateJobRunnerStatus(exitStatus: Int): Unit = {
    exitStatus match {
      case 0 => updateStatus(RunnerStatus.DONE)
      case _ => updateStatus(RunnerStatus.FAILED)
    }
  }

  /**
   * Attempts to get the status of a job by looking at if the finalExitStatus
   * promise has completed or not.
   * @return if the jobRunner has updated it's status or not.
   */
  def updateJobStatus(): Boolean = {
    if (finalExitStatus.isCompleted) {
      val completedExitStatus = finalExitStatus.future.value.get.get
      exitStatusUpdateJobRunnerStatus(completedExitStatus)
      true
    } else {
      // Make sure the status is update here, otherwise Queue will think
      // it's lots control over the job and kill it after 5 minutes.
      updateStatus(status)
      false
    }
  }
}
