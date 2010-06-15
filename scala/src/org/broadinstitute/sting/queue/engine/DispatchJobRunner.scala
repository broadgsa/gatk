package org.broadinstitute.sting.queue.engine

import edu.mit.broad.core.lsf.LocalLsfJob
import collection.JavaConversions._
import management.ManagementFactory
import java.io.File
import java.util.ArrayList
import org.broadinstitute.sting.queue.function.{DispatchFunction, QFunction}

trait DispatchJobRunner {
  type DispatchJobType
  private var dispatchJobs = Map.empty[DispatchFunction, DispatchJobType]

  protected def newJobName = DispatchJobRunner.nextJobName

  def dispatch(function: DispatchFunction, qGraph: QGraph)

  protected def addJob(function: DispatchFunction, dispatchJob: DispatchJobType) =
    dispatchJobs += function -> dispatchJob

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
        case qFunction: QFunction => previous :::= previousJobs(qFunction, qGraph)
      }
    }
    previous
  }
}

object DispatchJobRunner {
  private val jobNamePrefix = "Q-" + {
    var prefix = ManagementFactory.getRuntimeMXBean.getName
    val index = prefix.indexOf(".")
    if (index >= 0)
      prefix = prefix.substring(0, index)
    prefix
  }
  private var jobIndex = 0

  private def nextJobName = {
    jobIndex += 1
    jobNamePrefix + "-" + jobIndex
  }
}
