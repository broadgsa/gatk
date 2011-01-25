package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.ShellJob

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

    writeExec()
    job.shellScript = exec

    // Allow advanced users to update the job.
    updateJobRun(job)

    runStatus = RunnerStatus.RUNNING
    job.run()
    runStatus = RunnerStatus.DONE
  }

  def status = runStatus
}
