package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.{CommandLineFunction, QFunction}
import scala.collection.immutable.ListSet
import org.broadinstitute.sting.queue.util.IOUtils
import java.io.File

/**
 * Dispatches jobs to a compute cluster.
 */
trait DispatchJobRunner extends JobRunner {
  /** Type of the job. */
  type DispatchJobType
  /** An internal cache of all the jobs that have run by command line function. */
  private var dispatchJobs = Map.empty[CommandLineFunction, DispatchJobType]
  /** An internal list of functions that have no other dependencies. */
  private var waitJobsByGraph = Map.empty[QGraph, ListSet[DispatchJobType]]

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
  protected def previousJobs(function: CommandLineFunction, qGraph: QGraph) : List[DispatchJobType] =
    qGraph.previousJobs(function).map(dispatchJobs(_))

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
    var dirs = Set.empty[File]
    for (dir <- function.jobDirectories)
      dirs += IOUtils.dirLevel(dir, 2)
    if (dirs.size > 0)
      Some(dirs.mkString("cd ", " && cd ", ""))
    else
      None
  }
}
