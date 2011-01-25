package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Creates and stops Lsf706JobRunners
 */
class Lsf706JobManager extends JobManager[CommandLineFunction, Lsf706JobRunner] {
  def runnerType = classOf[Lsf706JobRunner]
  def functionType = classOf[CommandLineFunction]
  def create(function: CommandLineFunction) = new Lsf706JobRunner(function)
  override def tryStop(runners: List[Lsf706JobRunner]) { Lsf706JobRunner.tryStop(runners) }
}
