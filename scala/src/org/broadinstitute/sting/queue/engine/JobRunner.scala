package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Base interface for job runners.
 */
trait JobRunner {
  /**
   * Dispatches a function to the queue and returns immediately, unless the function is a DispatchWaitFunction
   * in which case it waits for all other terminal functions to complete.
   * @param function Command to run.
   * @param qGraph graph that holds the job, and if this is a dry run.
   */
  def run(function: CommandLineFunction, qGraph: QGraph)
}