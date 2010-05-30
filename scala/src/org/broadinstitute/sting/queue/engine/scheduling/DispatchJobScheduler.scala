package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.DirectedGraph
import edu.mit.broad.core.lsf.LocalLsfJob
import collection.JavaConversions._
import management.ManagementFactory
import java.io.File
import java.util.ArrayList
import org.broadinstitute.sting.queue.QArguments

/**
 * Dispatches jobs to LSF and then returns.
 */
class DispatchJobScheduler(jobGraph: DirectedGraph[ResourceNode, ResourceEdge], qArgs: QArguments)
        extends TopologicalJobScheduler(jobGraph, qArgs) {
  private var lsfJobs = Map.empty[ExecEdge, LocalLsfJob]
  private var lsfJobIndex = 0
  private val jvmName = ManagementFactory.getRuntimeMXBean.getName
  private val jobNamePrefix = "Q-" + jvmName

  def processExec(exec: ExecEdge) = {
    lsfJobIndex += 1
    val job = new LocalLsfJob
    val jobName = jobNamePrefix + "-" + lsfJobIndex
    val outputFile = jobName + ".out"
    val errorFile = jobName + ".err"
    val workingDir = lookup(exec, "jobWorkingDir", ".")
    val lsfProject = lookup(exec, "jobProject", "Queue")
    val queue = lookup(exec, "jobQueue", "broad")
    val memory = lookup(exec, "jobMemory", "2")

    var extraArgs = List("-r", "-R", "rusage[mem=" + memory + "]")

    val sourceJobs = sourceLsfJobs(exec)
    if (sourceJobs.size > 0) {
      extraArgs :::= List("-w", dependencyExpression(sourceJobs))
    }
    job.setName(jobName)
    job.setExtraBsubArgs(new ArrayList(extraArgs))
    job.setProject(lsfProject)
    job.setWorkingDir(new File(workingDir))
    job.setProject(lsfProject)
    job.setCommand(exec.commandString)
    job.setOutputFile(new File(workingDir, outputFile))
    job.setErrFile(new File(workingDir, errorFile))
    job.setQueue(queue)

    lsfJobs = lsfJobs.updated(exec, job)

    logger.info(job.getBsubCommand.mkString(" "))

    if (!qArgs.dryRun)
      job.start
  }

  /**
   * Walks up the graph looking for the previous LsfJobs for this node
   */
  private def sourceLsfJobs(edge: ResourceEdge) : List[LocalLsfJob] = {
    var sourceJobs = List.empty[LocalLsfJob]

    val source = this.jobGraph.getEdgeSource(edge)
    for (incomingEdge <- this.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

        // Stop recursing when we find a job along this edge and return it's job id
        case exec: ExecEdge => sourceJobs :::= List(lsfJobs(exec))

        // For any other type of edge find the LSF jobs preceeding the edge
        case resourceEdge: ResourceEdge => sourceJobs :::= sourceLsfJobs(resourceEdge)
      }
    }
    sourceJobs
  }

  private def dependencyExpression(jobs: List[LocalLsfJob]) = {
    jobs.toSet[LocalLsfJob].map(_.getName).mkString("ended(\"", "\") && ended(\"", "\")")
  }
}
