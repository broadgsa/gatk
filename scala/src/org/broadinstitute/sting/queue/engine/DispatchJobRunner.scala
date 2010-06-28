package org.broadinstitute.sting.queue.engine

import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.{DispatchFunction, QFunction}
import scala.collection.immutable.ListSet

trait DispatchJobRunner {
  type DispatchJobType
  private var dispatchJobs = Map.empty[DispatchFunction, DispatchJobType]
  private var waitJobsByGraph = Map.empty[QGraph, ListSet[DispatchJobType]]

  /**
   * Dispatches a function to the queue and returns immediately, unless the function is a DispatchWaitFunction
   * in which case it waits for all other terminal functions to complete.
   */
  def dispatch(function: DispatchFunction, qGraph: QGraph)

  protected def addJob(function: DispatchFunction, qGraph: QGraph,
      dispatchJob: DispatchJobType, previousJobs: List[DispatchJobType]) = {
    dispatchJobs += function -> dispatchJob
    var waitJobs = getWaitJobs(qGraph)
    for (previousJob <- previousJobs)
      waitJobs -= previousJob
    waitJobs += dispatchJob
    waitJobsByGraph += qGraph -> waitJobs
  }

  /**
   * Walks up the graph looking for the previous LsfJobs
   */
  protected def previousJobs(function: QFunction, qGraph: QGraph) : List[DispatchJobType] = {
    var previous = List.empty[DispatchJobType]

    val source = qGraph.jobGraph.getEdgeSource(function)
    for (incomingEdge <- qGraph.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

        // Stop recursing when we find a job along the edge and return its job id
        case dispatchFunction: DispatchFunction => previous :+= dispatchJobs(dispatchFunction)

        // For any other type of edge find the LSF jobs preceding the edge
        case qFunction: QFunction => previous = previousJobs(qFunction, qGraph) ::: previous
      }
    }
    previous
  }

  /**
   * Returns a set of jobs that have no following jobs in the graph.
   */
  protected def getWaitJobs(qGraph: QGraph) = {
    if (!waitJobsByGraph.contains(qGraph))
      waitJobsByGraph += qGraph -> ListSet.empty[DispatchJobType]
    waitJobsByGraph(qGraph)
  }
}
