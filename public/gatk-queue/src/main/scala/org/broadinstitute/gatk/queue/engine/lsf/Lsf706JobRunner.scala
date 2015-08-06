/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.engine.lsf

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.util._
import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.{LibLsf, LibBat}
import org.broadinstitute.gatk.utils.Utils
import org.broadinstitute.gatk.utils.jna.clibrary.LibC
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.LibBat.{submitReply, submit}
import org.broadinstitute.gatk.queue.engine.{RunnerStatus, CommandLineJobRunner}
import java.util.regex.Pattern
import java.lang.StringBuffer
import java.util.Date
import com.sun.jna.{Pointer, Structure, StringArray, NativeLong}
import com.sun.jna.ptr.IntByReference

/**
 * Runs jobs on an LSF compute cluster.
 */
class Lsf706JobRunner(val function: CommandLineFunction) extends CommandLineJobRunner with Logging {

  // Run the static initializer for Lsf706JobRunner
  Lsf706JobRunner

  /** Job Id of the currently executing job. */
  private var jobId = -1L
  override def jobIdString = jobId.toString

  protected override val minRunnerPriority = 1
  protected override val maxRunnerPriority = Lsf706JobRunner.maxUserPriority

  private val selectString = new StringBuffer()
  private val usageString = new StringBuffer()
  private val requestString = new StringBuffer()
  private val spanString = new StringBuffer()

  /**
   * Dispatches the function on the LSF cluster.
   */
  def start() {
    Lsf706JobRunner.lsfLibLock.synchronized {

      parseResourceRequest()

      val request = new submit
      for (i <- 0 until LibLsf.LSF_RLIM_NLIMITS)
        request.rLimits(i) = LibLsf.DEFAULT_RLIMIT;

      request.jobName = function.jobRunnerJobName.take(LibBat.MAX_JOB_NAME_LEN)
      request.options |= LibBat.SUB_JOB_NAME

      // Set the output file for stdout
      request.outFile = function.jobOutputFile.getPath
      request.options |= LibBat.SUB_OUT_FILE

      // Set the current working directory
      request.cwd = function.commandDirectory.getPath
      request.options3 |= LibBat.SUB3_CWD

      // If the error file is set specify the separate output for stderr
      if (function.jobErrorFile != null) {
        request.errFile = function.jobErrorFile.getPath
        request.options |= LibBat.SUB_ERR_FILE
      }

      // If the job queue is set specify the job queue
      if (function.jobQueue != null) {
        request.queue = function.jobQueue
        request.options |= LibBat.SUB_QUEUE
      }

      // If the resident set size is requested pass on the memory request
      if (function.residentRequest.isDefined) {
        val memInUnits = Lsf706JobRunner.convertUnits(function.residentRequest.get)
        appendRequest("select", selectString, "&&", "mem>%d".format(memInUnits))
        appendRequest("rusage", usageString, ",", "mem=%d".format(memInUnits))
      }

      //
      // Request multiple cores on the same host.  If nCoresRequest > 1, and we
      // aren't being jerks and stealing cores, set numProcessors and maxNumProcessors
      // and the span[host=1] parameters to get us exactly the right number of
      // cores on a single host
      //
      if ( function.nCoresRequest.getOrElse(1) > 1 ) {
        if ( function.qSettings.dontRequestMultipleCores )
          logger.warn("Sending multicore job %s to farm without requesting appropriate number of cores (%d)".format(
            function.shortDescription, function.nCoresRequest.get))
        else {
          request.numProcessors = function.nCoresRequest.get
          request.maxNumProcessors = request.numProcessors
          appendRequest("span", spanString, ",", "hosts=1")
        }
      }

      val resReq = getResourceRequest
      if (resReq.length > 0) {
        request.resReq = resReq
        request.options |= LibBat.SUB_RES_REQ
      }

      // If the resident set size limit is defined specify the memory limit
      if (function.residentLimit.isDefined) {
        val memInUnits = Lsf706JobRunner.convertUnits(function.residentLimit.get)
        request.rLimits(LibLsf.LSF_RLIMIT_RSS) = memInUnits
      }

      // If the priority is set (user specified Int) specify the priority
      val priority = functionPriority
      if (priority.isDefined) {
        request.userPriority = priority.get
        request.options2 |= LibBat.SUB2_JOB_PRIORITY
      }

      // Set the project to either the function or LSF default
      val project = if (function.jobProject != null) function.jobProject else Lsf706JobRunner.defaultProject
      if (project != null) {
        request.projectName = project
        request.options |= LibBat.SUB_PROJECT_NAME
      }

      // Set the esub names based on the job envorinment names
      if (!function.jobEnvironmentNames.isEmpty) {
        val argv = Array("", "-a", function.jobEnvironmentNames.mkString(" "))
        val setOptionResult = LibBat.setOption_(argv.length, new StringArray(argv), "a:", request, ~0, ~0, ~0, null);
        if (setOptionResult == -1)
          throw new QException("setOption_() returned -1 while setting esub");
      }

      if(!function.wallTime.isEmpty)
        request.rLimits(LibLsf.LSF_RLIMIT_RUN) = function.wallTime.get.toInt
      else
        // LSF specific: get the max runtime for the jobQueue and pass it for this job
        request.rLimits(LibLsf.LSF_RLIMIT_RUN) = Lsf706JobRunner.getRlimitRun(function.jobQueue)

      // Run the command as sh <jobScript>
      request.command = "sh " + jobScript

      // Allow advanced users to update the request via QFunction.updateJobRun()
      updateJobRun(request)

      updateStatus(RunnerStatus.RUNNING)
      Retry.attempt(() => {
        val reply = new submitReply
        jobId = LibBat.lsb_submit(request, reply)
        if (jobId < 0)
          throw new QException(LibBat.lsb_sperror("Unable to submit job"))
      }, 1, 5, 10)
      logger.info("Submitted LSF job id: " + jobId)
    }
  }

  override def checkUnknownStatus() {
    // TODO: Need a second pass through either of the two archive logs using lsb_geteventrecbyline() for disappeared jobs.
    // Can also tell if we wake up and the last time we saw status was greater than lsb_parameterinfo().cleanPeriod
    // LSB_SHAREDIR/cluster_name/logdir/lsb.acct (man bacct)
    // LSB_SHAREDIR/cluster_name/logdir/lsb.events (man bhist)
    logger.debug("Job Id %s status / exitStatus / exitInfo: ??? / ??? / ???".format(jobId))
    super.checkUnknownStatus()
  }

  private def parseResourceRequest() {
    requestString.setLength(0)
    selectString.setLength(0)
    usageString.setLength(0)
    spanString.setLength(0)

    requestString.append(function.jobResourceRequests.mkString(" "))
    extractSection(requestString, "select", selectString)
    extractSection(requestString, "rusage", usageString)
    extractSection(requestString, "span", spanString)
  }

  private def extractSection(requestString: StringBuffer, section: String, sectionString: StringBuffer) {
    val pattern = Pattern.compile(section + "\\s*\\[[^\\]]+\\]\\s*");
    val matcher = pattern.matcher(requestString.toString)
    if (matcher.find()) {
      sectionString.setLength(0)
      sectionString.append(matcher.group().trim())

      val sb = new StringBuffer
      matcher.appendReplacement(sb, "")
      matcher.appendTail(sb)

      requestString.setLength(0)
      requestString.append(sb)
    }
  }

  private def appendRequest(section: String, sectionString: StringBuffer, separator: String, request: String) {
    if (sectionString.length() == 0)
      sectionString.append(section).append("[").append(request).append("]")
    else
      sectionString.insert(sectionString.length() - 1, separator + request)
  }

  private def getResourceRequest = "%s %s %s %s".format(selectString, usageString, spanString, requestString).trim()
}

object Lsf706JobRunner extends Logging {
  private val lsfLibLock = new Object
  private val SIGTERM = 15

  /** Number of seconds for a non-normal exit status before we give up on expecting LSF to retry the function. */
  private val retryExpiredSeconds = 5 * 60

  /**
   * Initialize the Lsf library.
   */
  private val (defaultQueue, defaultProject, maxUserPriority) = {
    lsfLibLock.synchronized {
      if (LibBat.lsb_init("Queue") < 0)
        throw new QException(LibBat.lsb_sperror("lsb_init() failed"))

      val parameterInfo = LibBat.lsb_parameterinfo(null, null, 0);
      var defaultQueue: String = parameterInfo.defaultQueues
      val defaultProject = parameterInfo.defaultProject
      val maxUserPriority = parameterInfo.maxUserPriority

      if (defaultQueue != null && defaultQueue.indexOf(' ') > 0)
        defaultQueue = defaultQueue.split(" ")(0)

      (defaultQueue, defaultProject, maxUserPriority)
    }
  }

  /**
   * Bulk updates job statuses.
   * @param runners Runners to update.
   * @return runners which were updated.
   */
  def updateStatus(runners: Set[Lsf706JobRunner]) = {
    var updatedRunners = Set.empty[Lsf706JobRunner]

    Lsf706JobRunner.lsfLibLock.synchronized {
      val result = LibBat.lsb_openjobinfo(0L, null, null, null, null, LibBat.ALL_JOB)
      if (result < 0) {
        logger.error(LibBat.lsb_sperror("Unable to check LSF job info"))
      } else {
        try {
          val more = new IntByReference(result)
          while (more.getValue > 0) {
            val jobInfo = LibBat.lsb_readjobinfo(more)
            if (jobInfo == null) {
              logger.error(LibBat.lsb_sperror("Unable to read LSF job info"))
              more.setValue(0)
            } else {
              runners.find(runner => runner.jobId == jobInfo.jobId) match {
                case Some(runner) =>
                  updateRunnerStatus(runner, jobInfo)
                  updatedRunners += runner
                case None => /* not our job */
              }
            }
          }
        } finally {
          LibBat.lsb_closejobinfo()
        }
      }
    }

    updatedRunners
  }

  private def updateRunnerStatus(runner: Lsf706JobRunner, jobInfo: LibBat.jobInfoEnt) {
    val jobStatus = jobInfo.status
    val exitStatus = jobInfo.exitStatus
    val exitInfo = jobInfo.exitInfo
    val endTime = jobInfo.endTime

    logger.debug("Job Id %s status / exitStatus / exitInfo: 0x%02x / 0x%02x / 0x%02x".format(runner.jobId, jobStatus, exitStatus, exitInfo))

    def updateRunInfo() {
      // the platform LSF startTimes are in seconds, not milliseconds, so convert to the java convention
      runner.getRunInfo.startTime = new Date(jobInfo.startTime.longValue * 1000)
      runner.getRunInfo.doneTime = new Date(jobInfo.endTime.longValue * 1000)

      val exHostsList =
        if (jobInfo.numExHosts != 1) {
          // this is necessary because
          val exHostsString = "multipleHosts_" + jobInfo.numExHosts
          logger.debug("numExHosts = " + jobInfo.numExHosts + " != 1 for job " + runner.jobId + ", cannot safely get exhosts, setting to " + exHostsString)
          List(exHostsString)
        } else {
          jobInfo.exHosts.getStringArray(0).toSeq
        }

      //logger.warn("exHostsList = " + exHostsList)
      val exHosts = exHostsList.reduceLeft(_ + "," + _)
      //logger.warn("exHosts = " + exHosts)
      runner.getRunInfo.exechosts = exHosts
    }

    runner.updateStatus(
      if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_DONE)) {
        // Done successfully.
        updateRunInfo()
        RunnerStatus.DONE
      } else if (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_EXIT) && !willRetry(exitInfo, endTime)) {
        // Exited function that (probably) won't be retried.
        updateRunInfo()
        RunnerStatus.FAILED
      } else {
        // Note that we still saw the job in the system.
        RunnerStatus.RUNNING
      }
    )
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

  /**
   * Tries to stop any running jobs.
   * @param runners Runners to stop.
   */
  def tryStop(runners: Set[Lsf706JobRunner]) {
    lsfLibLock.synchronized {
      // lsb_killbulkjobs does not seem to forward SIGTERM,
      // only SIGKILL, so send the Ctrl-C (SIGTERM) one by one.
      for (runner <- runners.filterNot(_.jobId < 0)) {
        try {
          if (LibBat.lsb_signaljob(runner.jobId, SIGTERM) < 0)
            logger.error(LibBat.lsb_sperror("Unable to kill job " + runner.jobId))
        } catch {
          case e: Exception=>
            logger.error("Unable to kill job " + runner.jobId, e)
        }
      }
    }
  }

  /** The run limits for each queue. */
  private var queueRlimitRun = Map.empty[String,Int]

  /**
   * Returns the run limit in seconds for the queue.
   * If the queue name is null returns the length of the default queue.
   * @param queueName Name of the queue or null for the default queue.
   * @return the run limit in seconds for the queue.
   */
  private def getRlimitRun(queueName: String) = {
    lsfLibLock.synchronized {
      val queue = if (queueName == null) defaultQueue else queueName
      queueRlimitRun.get(queue) match {
        case Some(limit) => limit
        case None =>
          // Cache miss.  Go get the run limits from LSF.
          val queues = new StringArray(Array(queue))
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

  private lazy val unitDivisor: Double = {
    lsfLibLock.synchronized {
      val unitsParam: Array[LibLsf.config_param] = new LibLsf.config_param().toArray(2).asInstanceOf[Array[LibLsf.config_param]]
      unitsParam(0).paramName = "LSF_UNIT_FOR_LIMITS"

      Structure.autoWrite(unitsParam.asInstanceOf[Array[Structure]])
      if (LibLsf.ls_readconfenv(unitsParam(0), null) != 0)
        throw new QException(LibBat.lsb_sperror("ls_readconfenv() failed"))
      Structure.autoRead(unitsParam.asInstanceOf[Array[Structure]])

      unitsParam(0).paramValue match {
        case "MB" => 1 / 1024D
        case "GB" => 1D
        case "TB" => 1024D
        case "PB" => 1024D * 1024
        case "EB" => 1024D * 1024 * 1024
        case null => 1 / 1024D
      }
    }
  }

  private def convertUnits(gb: Double) = (gb / unitDivisor).ceil.toInt
}
