/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.queue.engine.slurm

import java.util.Date

import org.broadinstitute.gatk.queue.engine.parallelshell.ThreadSafeProcessController
import org.broadinstitute.gatk.queue.engine.{CommandLineJobRunner, RunnerStatus}
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.util.Logging
import org.broadinstitute.gatk.utils.Utils
import org.broadinstitute.gatk.utils.runtime.{OutputStreamSettings, ProcessSettings}

import scala.collection.mutable.ArrayBuffer
import scala.concurrent._
import scala.util.{Failure, Success}
import java.util.concurrent.Executors

/**
 * Runs multiple jobs locally without blocking.
 * Use this with care as it might not be the most efficient way to run things.
 * However, for some scenarios, such as running multiple single threaded
 * programs concurrently it can be quite useful.
 *
 * All this code is based on the normal shell runner in GATK Queue and all
 * credits for everything except the concurrency part goes to the GATK team.
 *
 * @author Peter van 't Hof
 *
 * @param function Command to run.
 */
class SlurmJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {

  // Controller on the thread that started the job
  val controller: ThreadSafeProcessController = new ThreadSafeProcessController()

  // Once the application exits this promise will be fulfilled.
  val finalExitStatus = Promise[Int]()

  /**
   * Runs the function on the local shell.
   */
  def start() {

    val srunCommand = new ArrayBuffer[String]
    srunCommand.append("srun")

    srunCommand.append(s"--job-name=${function.jobName}")
    function.qualityOfSerice.foreach(qos => srunCommand.append(s"--qos=$qos"))
    function.residentLimit.foreach(mem => srunCommand.append(s"--mem=${mem.ceil.toInt}G"))
    function.wallTime.foreach(time => srunCommand.append(s"--time=$time:00:00"))
    function.nCoresRequest.foreach(cores => srunCommand.append(s"--cpus-per-task=$cores"))
    function.jobNativeArgs.foreach(arg => srunCommand.append(arg.split(" "):_*))

    logger.info(s"Native arguments: ${srunCommand.tail.mkString(" ")}")

    jobScript.setExecutable(true)
    srunCommand.append(jobScript.getAbsolutePath)

    val commandLine = srunCommand.toArray
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

    implicit val ec =  SlurmJobRunner.ec

    // Run the command line process in a future.
    val executedFuture =
      future {
        this.function.waitBeforeJob.foreach(x => Thread.sleep(1000 * x))
        controller.exec(processSettings)
      }

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
      // Make sure the status is updated here, otherwise Queue will think
      // it has lost control over the job and kill it after 5 minutes.
      updateStatus(status)
      false
    }
  }
}

object SlurmJobRunner {
  // Spawn a new thread for each process
  val ec = new ExecutionContext {
    val threadPool = Executors.newFixedThreadPool(1000)

    def execute(runnable: Runnable) {
      threadPool.submit(runnable)
    }

    def reportFailure(t: Throwable) {}
  }
}