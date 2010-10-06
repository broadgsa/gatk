package org.broadinstitute.sting.queue.engine

/**
 * Base interface for job runners.
 */
trait JobRunner {
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
}
