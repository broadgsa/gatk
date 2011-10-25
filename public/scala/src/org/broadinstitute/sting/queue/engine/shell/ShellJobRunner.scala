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

package org.broadinstitute.sting.queue.engine.shell

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.engine.{RunnerStatus, CommandLineJobRunner}
import java.util.Date
import org.broadinstitute.sting.utils.Utils
import org.broadinstitute.sting.utils.runtime.{ProcessSettings, OutputStreamSettings, ProcessController}

/**
 * Runs jobs one at a time locally
 */
class ShellJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner {
  // Controller on the thread that started the job
  private var controller: ProcessController = null

  /**
   * Runs the function on the local shell.
   * @param function Command to run.
   */
  def start() {
    val commandLine = Array("sh", jobScript.getAbsolutePath)
    val stdoutSettings = new OutputStreamSettings
    val stderrSettings = new OutputStreamSettings
    val mergeError = (function.jobErrorFile != null)

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
    controller = ProcessController.getThreadLocal
    val exitStatus = controller.exec(processSettings).getExitValue
    getRunInfo.doneTime = new Date()
    updateStatus(if (exitStatus == 0) RunnerStatus.DONE else RunnerStatus.FAILED)
  }

  /**
   * Possibly invoked from a shutdown thread, find and
   * stop the controller from the originating thread
   */
  def tryStop() {
    // Assumes that after being set the job may be
    // reassigned but will not be reset back to null
    if (controller != null) {
      try {
        controller.tryDestroy()
      } catch {
        case e =>
          logger.error("Unable to kill shell job: " + function.description)
      }
    }
  }
}
