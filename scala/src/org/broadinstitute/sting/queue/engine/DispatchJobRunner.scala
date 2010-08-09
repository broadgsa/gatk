package org.broadinstitute.sting.queue.engine

import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.{CommandLineFunction, QFunction}
import scala.collection.immutable.ListSet

/**
 * Dispatches jobs to a compute cluster.
 */
trait DispatchJobRunner {
  /** Type of the job. */
  type DispatchJobType
  /** An internal cache of all the jobs that have run by command line function. */
  private var dispatchJobs = Map.empty[CommandLineFunction, DispatchJobType]
  /** An internal list of functions that have no other dependencies. */
  private var waitJobsByGraph = Map.empty[QGraph, ListSet[DispatchJobType]]

  /**
   * Dispatches a function to the queue and returns immediately, unless the function is a DispatchWaitFunction
   * in which case it waits for all other terminal functions to complete.
   * @param function Command to run.
   * @param qGraph graph that holds the job, and if this is a dry run.
   */
  def dispatch(function: CommandLineFunction, qGraph: QGraph)

  /**
   * Adds the job to the internal cache of previous jobs and removes the previous jobs that
   * the job was dependent on from the list of function that have no dependencies.
   * @param function CommandLineFunction to add to the list.
   * @param qGraph Current qGraph being iterated over.
   * @param dispatchJob The job that is being added to the cache.
   * @param previousJobs The previous jobs that the job was dependent one.
   */
  protected def addJob(function: CommandLineFunction, qGraph: QGraph,
      dispatchJob: DispatchJobType, previousJobs: Iterable[DispatchJobType]) = {
    dispatchJobs += function -> dispatchJob
    var waitJobs = getWaitJobs(qGraph)
    for (previousJob <- previousJobs)
      waitJobs -= previousJob
    waitJobs += dispatchJob
    waitJobsByGraph += qGraph -> waitJobs
  }

  /**
   * Walks up the graph looking for the previous LsfJobs.
   * @param function Function to examine for a previous command line job.
   * @param qGraph The graph that contains the jobs.
   * @return A list of prior jobs.
   */
  protected def previousJobs(function: QFunction, qGraph: QGraph) : List[DispatchJobType] = {
    var previous = List.empty[DispatchJobType]

    val source = qGraph.jobGraph.getEdgeSource(function)
    for (incomingEdge <- qGraph.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

        // Stop recursing when we find a job along the edge and return its job id
        case dispatchFunction: CommandLineFunction => previous :+= dispatchJobs(dispatchFunction)

        // For any other type of edge find the LSF jobs preceding the edge
        case qFunction: QFunction => previous ++= previousJobs(qFunction, qGraph)
      }
    }
    previous
  }

  /**
   * Returns a set of jobs that have no following jobs in the graph.
   * @param qGraph The graph that contains the jobs.
   * @return ListSet[DispatchJobType] of previous jobs that have no dependent jobs.
   */
  protected def getWaitJobs(qGraph: QGraph) = {
    if (!waitJobsByGraph.contains(qGraph))
      waitJobsByGraph += qGraph -> ListSet.empty[DispatchJobType]
    waitJobsByGraph(qGraph)
  }

  /**
   * Builds a command line that can be run to force an automount of the directories.
   * @param function Function to look jobDirectories.
   * @return A "cd <dir_1> [&& cd <dir_n>]" command.
   */
  protected def mountCommand(function: CommandLineFunction) = {
    val dirs = function.jobDirectories
    if (dirs.size > 0)
      Some("\'" + dirs.mkString("cd ", " && cd ", "") + "\'")
    else
      None
  }
}
