package org.broadinstitute.sting.queue.engine.scheduling

import org.jgrapht.graph.SimpleDirectedGraph
import edu.mit.broad.core.lsf.LocalLsfJob
import collection.JavaConversions._
import management.ManagementFactory
import org.broadinstitute.sting.queue.QException
import java.io.File
import java.util.{ArrayList, Properties}
import org.broadinstitute.sting.queue.engine.Pipeline

/**
 * Dispatches jobs to LSF and then returns.
 */
class DispatchJobScheduler(jobGraph: SimpleDirectedGraph[ResourceNode, ResourceEdge]) extends TopologicalJobScheduler(jobGraph) {
  private var lsfJobs = Map.empty[ExecEdge, LocalLsfJob]
  private var lsfJobIndex = 0
  private val jvmName = ManagementFactory.getRuntimeMXBean.getName
  private val jobNamePrefix = "Q-" + jvmName

  protected def traversedExec(exec: ExecEdge) = {
    lsfJobIndex += 1
    val props = new Properties
    exec.args.foreach(x => props.put(x._1, x._2))
    val job = new LocalLsfJob
    val jobName = jobNamePrefix + "-" + lsfJobIndex
    val outputFile = jobName + ".out"
    val errorFile = jobName + ".err"
    val workingDir = props.getProperty("jobWorkingDir", ".")
    val lsfProject = props.getProperty("jobProject", "Queue")
    val queue = props.getProperty("jobQueue", "broad")
    val memory = props.getProperty("jobMemory", "2")

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

    if (!Pipeline.dryRun)
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

        // Throw error for a new edge type that we don't know how to handle
        case default => throw new QException("Unknown edge type: " + default)
      }
    }
    sourceJobs
  }

  private def dependencyExpression(jobs: List[LocalLsfJob]) = {
    jobs.toSet[LocalLsfJob].map(_.getName).mkString("done(\"", "\") && done(\"", "\")")
  }
}
