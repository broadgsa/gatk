package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base interface for job runners.
 */
trait JobRunner[TFunction <: QFunction] {
  /**
   * Runs the function.
   * After the function returns the status of the function should
   * be RUNNING, FAILED, or DONE.
   * @param function Command to run.
   */
  def start()

  /**
   * Returns the current run status.
   * Must only be called AFTER start().
   * @return RUNNING, DONE, or FAILED.
   */
  def status: RunnerStatus.Value

  /**
   * Returns the function to be run.
   */
  def function: TFunction

  /**
   * Removes all temporary files used for this job.
   */
  def removeTemporaryFiles() {
  }

  /**
   * Calls back to a hook that an expert user can setup to modify a job.
   * @param value Value to modify.
   */
  protected def updateJobRun(value: Any) {
    val updater = function.updateJobRun
    if (updater != null)
      if (updater.isDefinedAt(value))
        updater(value)
  }
}
