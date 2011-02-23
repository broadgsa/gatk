package org.broadinstitute.sting.queue.engine

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.jna.lsf.v7_0_6.{LibLsf, LibBat}
import org.broadinstitute.sting.utils.Utils
import org.broadinstitute.sting.jna.clibrary.LibC
import java.util.Date
import org.broadinstitute.sting.jna.lsf.v7_0_6.LibBat.{submitReply, submit}
import com.sun.jna.ptr.IntByReference
import com.sun.jna.{StringArray, NativeLong}

/**
 * Runs jobs on an LSF compute cluster.
 */
class Lsf706JobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {

  // Run the static initializer for Lsf706JobRunner
  Lsf706JobRunner

  /** Job Id of the currently executing job. */
  var jobId = -1L

  /** Last known run status */
  private var runStatus: RunnerStatus.Value = _

  /**
   * Dispatches the function on the LSF cluster.
   * @param function Command to run.
   */
  def start() = {
    Lsf706JobRunner.lsfLibLock.synchronized {
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

      if (function.jobPriority.isDefined) {
        request.userPriority = function.jobPriority.get
        request.options2 |= LibBat.SUB2_JOB_PRIORITY
      }

      request.rLimits(LibLsf.LSF_RLIMIT_RUN) = Lsf706JobRunner.getRlimitRun(function.jobQueue)

      writeExec()
      request.command = "sh " + exec

      // Allow advanced users to update the request.
      updateJobRun(request)

      runStatus = RunnerStatus.RUNNING
      Retry.attempt(() => {
        val reply = new submitReply
        jobId = LibBat.lsb_submit(request, reply)
        if (jobId < 0)
          throw new QException(LibBat.lsb_sperror("Unable to submit job"))
      }, 1, 5, 10)
      logger.info("Submitted LSF job id: " + jobId)
    }
  }

  /**
   * Updates and returns the status.
   */
  def status = {
    Lsf706JobRunner.lsfLibLock.synchronized {
      var jobStatus = LibBat.JOB_STAT_UNKWN
      var exitStatus = 0
      var exitInfo = 0
      var endTime: NativeLong = null

      var result = 0
      Retry.attempt(() => {
        result = LibBat.lsb_openjobinfo(jobId, null, null, null, null, LibBat.ALL_JOB)
        if (result < 0)
          throw new QException(LibBat.lsb_sperror("Unable to open LSF job info for job id: " + jobId))
      }, 0.5, 1, 2)
      try {
        if (result > 1)
          throw new QException(LibBat.lsb_sperror("Recieved " + result + " LSF results for job id: " + jobId))
        else if (result == 1) {
          val more = new IntByReference(result)
          val jobInfo = LibBat.lsb_readjobinfo(more)
          if (jobInfo == null)
            throw new QException(LibBat.lsb_sperror("lsb_readjobinfo returned null for job id: " + jobId))
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
          logger.error("Unknown status for %d seconds: job id %d: %s".format(unknownStatusMaxSeconds, jobId, function.description))
        }
      } else {
        // Reset the last time an unknown status was seen.
        firstUnknownTime = None

        if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_EXIT) && !willRetry(exitInfo, endTime)) {
          // Exited function that (probably) won't be retried.
          runStatus = RunnerStatus.FAILED
        } else if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_DONE)) {
          // Done successfully.
          runStatus = RunnerStatus.DONE
        }
      }

      runStatus
    }
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
  private val lsfLibLock = new Object
  private val SIGTERM = 15

  init()

  /** The name of the default queue. */
  private var defaultQueue: String = _

  /** The run limits for each queue. */
  private var queueRlimitRun = Map.empty[String,Int]

  /**
   * Initialize the Lsf library.
   */
  private def init() = {
    lsfLibLock.synchronized {
      if (LibBat.lsb_init("Queue") < 0)
        throw new QException(LibBat.lsb_sperror("lsb_init() failed"))
    }
  }

  /**
   * Returns the run limit in seconds for the queue.
   * If the queue name is null returns the length of the default queue.
   * @param queue Name of the queue or null for the default queue.
   * @return the run limit in seconds for the queue.
   */
  def getRlimitRun(queue: String) = {
    lsfLibLock.synchronized {
      if (queue == null) {
        if (defaultQueue != null) {
          queueRlimitRun(defaultQueue)
        } else {
          // Get the info on the default queue.
          val numQueues = new IntByReference(1)
          val queueInfo = LibBat.lsb_queueinfo(null, numQueues, null, null, 0)
          if (queueInfo == null)
            throw new QException(LibBat.lsb_sperror("Unable to get LSF queue info for the default queue"))
          defaultQueue = queueInfo.queue
          val limit = queueInfo.rLimits(LibLsf.LSF_RLIMIT_RUN)
          queueRlimitRun += defaultQueue -> limit
          limit
        }
      } else {
        queueRlimitRun.get(queue) match {
          case Some(limit) => limit
          case None =>
          // Cache miss.  Go get the run limits from LSF.
            val queues = new StringArray(Array[String](queue))
            val numQueues = new IntByReference(1)
            val queueInfo = LibBat.lsb_queueinfo(queues, numQueues, null, null, 0)
            if (queueInfo == null)
              throw new QException(LibBat.lsb_sperror("Unable to get LSF queue info for queue: " + queue))
            val limit = queueInfo.rLimits(LibLsf.LSF_RLIMIT_RUN)
            queueRlimitRun += queue -> limit
            limit
        }
      }
    }
  }

  /**
   * Tries to stop any running jobs.
   * @param runners Runners to stop.
   */
  def tryStop(runners: List[Lsf706JobRunner]) {
    lsfLibLock.synchronized {
      // lsb_killbulkjobs does not seem to forward SIGTERM,
      // only SIGKILL, so send the Ctrl-C (SIGTERM) one by one.
      for (jobRunner <- runners.filterNot(_.jobId < 0)) {
        try {
          if (LibBat.lsb_signaljob(jobRunner.jobId, SIGTERM) < 0)
            logger.error(LibBat.lsb_sperror("Unable to kill job " + jobRunner.jobId))
        } catch {
          case e =>
            logger.error("Unable to kill job " + jobRunner.jobId, e)
        }
      }
    }
  }
}
