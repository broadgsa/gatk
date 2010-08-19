package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.{CommandLineFunction, DispatchWaitFunction}
import org.broadinstitute.sting.queue.util.{IOUtils, LsfJob, Logging}

/**
 * Runs jobs on an LSF compute cluster.
 */
class LsfJobRunner extends DispatchJobRunner with Logging {
  type DispatchJobType = LsfJob

  /**
   * Dispatches the function on the LSF cluster.
   * @param function Command to run.
   * @param qGraph graph that holds the job, and if this is a dry run.
   */
  def run(function: CommandLineFunction, qGraph: QGraph) = {
    val job = new LsfJob
    job.name = function.jobName
    job.outputFile = function.jobOutputFile
    job.errorFile = function.jobErrorFile
    job.project = function.jobProject
    job.queue = function.jobQueue
    job.command = function.commandLine

    if (!IOUtils.CURRENT_DIR.getCanonicalFile.equals(function.commandDirectory))
      job.workingDir = function.commandDirectory

    if (function.jobRestartable)
      job.extraBsubArgs :+= "-r"

    if (function.memoryLimit.isDefined)
      job.extraBsubArgs ++= List("-R", "rusage[mem=" + function.memoryLimit.get + "]")

    val previous: Iterable[LsfJob] =
      if (function.isInstanceOf[DispatchWaitFunction]) {
        job.waitForCompletion = true
        getWaitJobs(qGraph)
      } else {
        previousJobs(function, qGraph)
      }

    mountCommand(function) match {
      case Some(command) => job.preExecCommand = command
      case None => /* ignore */
    }

    if (previous.size > 0)
      job.extraBsubArgs ++= List("-w", dependencyExpression(previous, function.jobRunOnlyIfPreviousSucceed))

    addJob(function, qGraph, job, previous)

    if (logger.isDebugEnabled) {
      logger.debug(function.commandDirectory + " > " + job.bsubCommand.mkString(" "))
    } else {
      logger.info(job.bsubCommand.mkString(" "))
    }

    if (!qGraph.dryRun)
      job.run
  }

  /**
   * Returns the dependency expression for the prior jobs.
   * @param jobs Previous jobs this job is dependent on.
   * @param runOnSuccess Run the job only if the previous jobs succeed.
   * @return The dependency expression for the prior jobs.
   */
  private def dependencyExpression(jobs: Iterable[LsfJob], runOnSuccess: Boolean) = {
    val jobIds = jobs.toSet[LsfJob].map(_.bsubJobId)
    if (runOnSuccess)
      jobIds.mkString("done(", ") && done(", ")")
    else
      jobIds.mkString("ended(", ") && ended(", ")")
  }
}
