package org.broadinstitute.sting.queue.engine

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.jna.lsf.v7_0_6.{LibLsf, LibBat}
import org.broadinstitute.sting.jna.lsf.v7_0_6.LibBat.{signalBulkJobs, submitReply, submit}
import org.broadinstitute.sting.utils.Utils
import com.sun.jna.{NativeLong, Memory}
import org.broadinstitute.sting.jna.clibrary.LibC
import java.util.Date

/**
 * Runs jobs on an LSF compute cluster.
 */
class Lsf706JobRunner(function: CommandLineFunction) extends LsfJobRunner(function) with Logging {

  // Run the static initializer for Lsf706JobRunner
  Lsf706JobRunner

  /**
   * Dispatches the function on the LSF cluster.
   * @param function Command to run.
   */
  def start() = {
    try {
      val request = new submit
      for (i <- 0 until LibLsf.LSF_RLIM_NLIMITS)
          request.rLimits(i) = LibLsf.DEFAULT_RLIMIT;

      request.outFile = function.jobOutputFile.getPath
      request.options |= LibBat.SUB_OUT_FILE

      if (function.jobErrorFile != null) {
        request.errFile = function.jobErrorFile.getPath
        request.options |= LibBat.SUB_ERR_FILE
      }

      if (function.jobProject != null) {
        request.projectName = function.jobProject
        request.options |= LibBat.SUB_PROJECT_NAME
      }

      if (function.jobQueue != null) {
        request.queue = function.jobQueue
        request.options |= LibBat.SUB_QUEUE
      }

      if (IOUtils.absolute(new File(".")) != function.commandDirectory) {
        request.cwd = function.commandDirectory.getPath
        request.options3 |= LibBat.SUB3_CWD
      }

      if (function.jobRestartable) {
        request.options |= LibBat.SUB_RERUNNABLE
      }

      if (function.memoryLimit.isDefined) {
        request.resReq = "rusage[mem=" + function.memoryLimit.get + "]"
        request.options |= LibBat.SUB_RES_REQ
      }

      if (function.description != null) {
        request.jobName = function.description.take(1000)
        request.options |= LibBat.SUB_JOB_NAME
      }

      if (function.jobLimitSeconds.isDefined) {
        request.rLimits(LibLsf.LSF_RLIMIT_RUN) = function.jobLimitSeconds.get
      }

      writeExec()
      request.command = "sh " + exec

      // Allow advanced users to update the request.
      updateJobRun(request)

      if (logger.isDebugEnabled) {
        logger.debug("Starting: " + function.commandDirectory + " > " + bsubCommand)
      } else {
        logger.info("Starting: " + bsubCommand)
      }

      function.deleteLogs()
      function.deleteOutputs()
      function.mkOutputDirectories()

      runStatus = RunnerStatus.RUNNING
      Retry.attempt(() => {
        val reply = new submitReply
        jobId = LibBat.lsb_submit(request, reply)
        if (jobId < 0)
          throw new QException(LibBat.lsb_sperror("Unable to submit job"))
      }, 1, 5, 10)
      logger.info("Submitted LSF job id: " + jobId)
    } catch {
      case e =>
        runStatus = RunnerStatus.FAILED
        try {
          removeTemporaryFiles()
          function.failOutputs.foreach(_.createNewFile())
          writeStackTrace(e)
        } catch {
          case _ => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + bsubCommand, e)
    }
  }

  /**
   * Updates and returns the status.
   */
  def status = {
    try {
      var jobStatus = LibBat.JOB_STAT_NULL
      var exitStatus = 0
      var exitInfo = 0
      var endTime: NativeLong = null

      LibBat.lsb_openjobinfo(jobId, null, null, null, null, LibBat.ALL_JOB)
      try {
        val jobInfo = LibBat.lsb_readjobinfo(null)
        if (jobInfo == null) {
          jobStatus = LibBat.JOB_STAT_UNKWN
          exitStatus = 0
          exitInfo = 0
          endTime = null
        } else {
          jobStatus = jobInfo.status
          exitStatus = jobInfo.exitStatus
          exitInfo = jobInfo.exitInfo
          endTime = jobInfo.endTime
        }
      } finally {
        LibBat.lsb_closejobinfo()
      }

      logger.debug("Job Id %s status / exitStatus / exitInfo: 0x%02x / 0x%02x / 0x%02x".format(jobId, jobStatus, exitStatus, exitInfo))

      if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_UNKWN)) {
        val now = new Date().getTime

        if (firstUnknownTime.isEmpty) {
          firstUnknownTime = Some(now)
          logger.debug("First unknown status for job id: " + jobId)
        }

        if ((firstUnknownTime.get - now) >= (unknownStatusMaxSeconds * 1000L)) {
          // Unknown status has been returned for a while now.
          runStatus = RunnerStatus.FAILED
          try {
            removeTemporaryFiles()
            function.failOutputs.foreach(_.createNewFile())
          } catch {
            case _ => /* ignore errors in the error handler */
          }
          logger.error("Error: " + bsubCommand + ", unknown status for " + unknownStatusMaxSeconds + " seconds.")
        }
      } else {
        // Reset the last time an unknown status was seen.
        firstUnknownTime = None

        if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_EXIT) && !willRetry(exitInfo, endTime)) {
          // Exited function that (probably) won't be retried.
          runStatus = RunnerStatus.FAILED
          try {
            removeTemporaryFiles()
            function.failOutputs.foreach(_.createNewFile())
          } catch {
            case _ => /* ignore errors in the error handler */
          }
          logger.error("Error: " + bsubCommand)
          tailError()
        } else if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_DONE)) {
          // Done successfully.
          removeTemporaryFiles()
          function.doneOutputs.foreach(_.createNewFile())
          runStatus = RunnerStatus.DONE
          logger.info("Done: " + bsubCommand)
        }
      }
    } catch {
      case e =>
        runStatus = RunnerStatus.FAILED
        try {
          removeTemporaryFiles()
          function.failOutputs.foreach(_.createNewFile())
          writeStackTrace(e)
        } catch {
          case _ => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + bsubCommand, e)
    }

    runStatus
  }

  /**
   * Returns true if LSF is expected to retry running the function.
   * @param exitInfo The reason the job exited.
   * @param endTime THe time the job exited.
   * @return true if LSF is expected to retry running the function.
   */
  private def willRetry(exitInfo: Int, endTime: NativeLong) = {
    exitInfo match {
      case LibBat.EXIT_NORMAL => false
      case _ => {
        val seconds = LibC.difftime(LibC.time(null), endTime)
        (seconds <= retryExpiredSeconds)
      }
    }
  }

}

object Lsf706JobRunner extends Logging {
  init()

  /**
   * Initialize the Lsf library.
   */
  private def init() = {
    if (LibBat.lsb_init("Queue") < 0)
      throw new QException(LibBat.lsb_sperror("lsb_init() failed"))
  }

  /**
   * Tries to stop any running jobs.
   * @param runners Runners to stop.
   */
  def tryStop(runners: List[JobRunner[_]]) = {
    val lsfJobRunners = runners.filter(_.isInstanceOf[Lsf706JobRunner]).map(_.asInstanceOf[Lsf706JobRunner])
    if (lsfJobRunners.size > 0) {
      for (jobRunners <- lsfJobRunners.filterNot(_.jobId < 0).grouped(10)) {
        try {
          val njobs = jobRunners.size
          val signalJobs = new signalBulkJobs
          signalJobs.jobs = {
            val p = new Memory(8 * njobs)
            p.write(0, jobRunners.map(_.jobId).toArray, 0, njobs)
            p
          }
          signalJobs.njobs = njobs
          signalJobs.signal = 9

          if (LibBat.lsb_killbulkjobs(signalJobs) < 0)
            throw new QException(LibBat.lsb_sperror("lsb_killbulkjobs failed"))
        } catch {
          case e =>
            logger.error("Unable to kill all jobs.", e)
        }
        try {
          jobRunners.foreach(_.removeTemporaryFiles())
        } catch {
          case e => /* ignore */
        }
      }
    }
  }
}
