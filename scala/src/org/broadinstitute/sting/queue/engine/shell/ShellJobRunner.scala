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
import org.broadinstitute.sting.queue.util.ShellJob
import org.broadinstitute.sting.queue.engine.{RunnerStatus, CommandLineJobRunner}

/**
 * Runs jobs one at a time locally
 */
class ShellJobRunner(val function: CommandLineFunction) extends CommandLineJobRunner {
  private var runStatus: RunnerStatus.Value = _

  /**
   * Runs the function on the local shell.
   * @param function Command to run.
   */
  def start() {
    val job = new ShellJob

    job.workingDir = function.commandDirectory
    job.outputFile = function.jobOutputFile
    job.errorFile = function.jobErrorFile

    job.shellScript = jobScript

    // Allow advanced users to update the job.
    updateJobRun(job)

    runStatus = RunnerStatus.RUNNING
    job.run()
    runStatus = RunnerStatus.DONE
  }

  def status = runStatus
}
