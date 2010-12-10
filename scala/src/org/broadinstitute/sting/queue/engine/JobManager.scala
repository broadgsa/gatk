package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.QFunction

/**
 * Creates and stops JobRunners
 */
trait JobManager[TFunction <: QFunction, TRunner <: JobRunner[TFunction]] {
  def create(function: TFunction): TRunner
  def tryStop(runners: List[JobRunner[_]]) = {}
}
