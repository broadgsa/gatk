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

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.engine.{RunnerStatus, CommandLineJobRunner}

class GridEngineJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {
  // Run the static initializer for GridEngineJobRunner
  GridEngineJobRunner

  def start() = {
    // TODO: Copy settings from function to GridEngine syntax.
    /*
    val gridEngineJob = new ...

    // Set the display name to 4000 characters of the description (or whatever the GE max is)
    gridEngineJob.displayName = function.description.take(4000)

    // Set the output file for stdout
    gridEngineJob.outputFile = function.jobOutputFile.getPath

    // Set the current working directory
    gridEngineJob.workingDirectory = function.commandDirectory.getPath

    // If the error file is set specify the separate output for stderr
    if (function.jobErrorFile != null) {
      gridEngineJob.errFile = function.jobErrorFile.getPath
    }

    // If a project name is set specify the project name
    if (function.jobProject != null) {
      gridEngineJob.projectName = function.jobProject
    }

    // If the job queue is set specify the job queue
    if (function.jobQueue != null) {
      gridEngineJob.queue = function.jobQueue
    }

    // If the memory limit is set (GB) specify the memory limit
    if (function.memoryLimit.isDefined) {
      gridEngineJob.jobMemoryLimit = function.memoryLimit.get + "GB"
    }

    // If the priority is set (user specified Int) specify the priority
    if (function.jobPriority.isDefined) {
      gridEngineJob.jobPriority = function.jobPriority.get
    }

    // Instead of running the function.commandLine, run "sh <jobScript>"
    gridEngineJob.command = "sh " + jobScript

    // Store the status so it can be returned in the status method.
    myStatus = RunnerStatus.RUNNING

    // Start the job and store the id so it can be killed in tryStop
    myJobId = gridEngineJob.start()
    */

    logger.warn("TODO: implement Grid Engine support")
  }

  // TODO: Return the latest status: RUNNING, FAILED, or DONE
  def status = throw new RuntimeException("TODO: Grid Engine return status such as: " + RunnerStatus.FAILED)
}

object GridEngineJobRunner extends Logging {
  initGridEngine()

  /**
   * Initialize the Grid Engine library.
   */
  private def initGridEngine() {
    // TODO: Init
    logger.warn("TODO Grid Engine: Initialize here.")
  }

  /**
   * Updates the status of a list of jobs.
   * @param runners Runners to update.
   */
  def updateStatus(runners: Set[GridEngineJobRunner]) {
    // TODO: Bulk update. If not possible this method can be removed here and in GridEngineJobManager.
  }

  /**
   * Tries to stop any running jobs.
   * @param runners Runners to stop.
   */
  def tryStop(runners: Set[GridEngineJobRunner]) {
    // TODO: Stop runners. SIGTERM(15) is preferred to SIGKILL(9).
  }
}
