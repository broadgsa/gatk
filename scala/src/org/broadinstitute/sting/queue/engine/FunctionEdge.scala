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

  def status = {
    if (currentStatus == RunnerStatus.PENDING || currentStatus == RunnerStatus.RUNNING)
      if (runner != null)
        currentStatus = runner.status
    currentStatus
  }

  def resetToPending() = {
    currentStatus = RunnerStatus.PENDING
    function.doneOutputs.foreach(_.delete())
    function.failOutputs.foreach(_.delete())
  }

  def inputs = function.inputs
  def outputs = function.outputs
  override def dotString = function.dotString
}
