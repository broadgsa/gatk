package org.broadinstitute.sting.queue.engine

import collection.JavaConversions._
import edu.mit.broad.core.lsf.LocalLsfJob
import org.broadinstitute.sting.queue.function.DispatchFunction
import java.util.ArrayList
import org.broadinstitute.sting.queue.util.Logging

trait LsfJobRunner extends DispatchJobRunner with Logging {
  type DispatchJobType = LocalLsfJob

  def dispatch(function: DispatchFunction, qGraph: QGraph) = {
    val job = new LocalLsfJob
    job.setName(function.jobName)
    job.setOutputFile(function.jobOutputFile)
    job.setErrFile(function.jobErrorFile)
    job.setWorkingDir(function.commandDirectory)
    job.setProject(function.jobProject)
    job.setQueue(function.jobQueue)
    job.setCommand(function.commandLine)

    var extraArgs = List("-r")

    if (function.memoryLimit.isDefined)
      extraArgs :::= List("-R", "rusage[mem=" + function.memoryLimit.get + "]")

    val previous = previousJobs(function, qGraph)
    if (previous.size > 0)
      extraArgs :::= List("-w", dependencyExpression(previous))

    job.setExtraBsubArgs(new ArrayList(extraArgs))

    addJob(function, job)

    if (logger.isDebugEnabled) {
      logger.debug(function.commandDirectory + " > " + job.getBsubCommand.mkString(" "))
    } else {
      logger.info(job.getBsubCommand.mkString(" "))
    }

    if (!qGraph.dryRun)
      job.start
  }

  private def dependencyExpression(jobs: List[LocalLsfJob]) = {
    jobs.toSet[LocalLsfJob].map(_.getName).mkString("ended(\"", "\") && ended(\"", "\")")
  }
}
