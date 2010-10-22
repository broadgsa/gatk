package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.QFunction

/**
 * An edge in the QGraph that runs a QFunction.
 * The edge is created first to determine inter-node dependencies,
 * and then the runner is specified later when the time comes to
 * execute the function in the edge.
 */
class FunctionEdge(var function: QFunction) extends QEdge {
  var runner: JobRunner =_

  /**
   * The number of times this edge has been run.
   */
  var retries = 0

  /**
   * Initializes with the current status of the function.
   */
  private var currentStatus = {
    val isDone = function.isDone
    val isFail = function.isFail
    if (isFail.isDefined && isFail.get)
      RunnerStatus.FAILED
    else if (isDone.isDefined && isDone.get)
      RunnerStatus.DONE
    else
      RunnerStatus.PENDING
  }

  /**
   * Returns the current status of the edge.
   */
  def status = {
    if (currentStatus == RunnerStatus.PENDING || currentStatus == RunnerStatus.RUNNING)
      if (runner != null)
        currentStatus = runner.status
    currentStatus
  }

  /**
   * Marks this edge as skipped as it is not needed for the current run.
   */
  def markAsSkipped() = {
    currentStatus = RunnerStatus.SKIPPED
  }

  /**
   * Resets the edge to pending status.
   */
  def resetToPending(cleanOutputs: Boolean) = {
    currentStatus = RunnerStatus.PENDING
    if (cleanOutputs)
      function.deleteOutputs()
    runner = null
  }

  def inputs = function.inputs
  def outputs = function.outputs
  override def dotString = function.dotString
}
