package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Creates and stops Lsf706JobRunners
 */
class Lsf706JobManager extends JobManager[CommandLineFunction, Lsf706JobRunner] {
  def create(function: CommandLineFunction) = new Lsf706JobRunner(function)
  override def tryStop(runners: List[JobRunner[_]]) = Lsf706JobRunner.tryStop(runners)
}
