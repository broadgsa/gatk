package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.QFunction

/**
 * Creates and stops JobRunners
 */
trait JobManager[TFunction <: QFunction, TRunner <: JobRunner[TFunction]] {
  /** The class type of the runner.  Available at runtime even after erasure. */
  def functionType: Class[TFunction]

  /** The class type of the functions processed by the runner.  Available at runtime even after erasure. */
  def runnerType: Class[TRunner]

  /** Creates a new runner.
   * @param function Function for the runner.
   */
  def create(function: TFunction): TRunner

  /**
   * Stops a list of functions.
   * @param runner Runners to stop.
   */
  def tryStop(runners: List[TRunner]) {
  }
}
