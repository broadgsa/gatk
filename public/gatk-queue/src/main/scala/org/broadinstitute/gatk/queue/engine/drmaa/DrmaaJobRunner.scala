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

package org.broadinstitute.gatk.queue.engine.drmaa

import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.queue.util.{Logging,Retry}
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.{RunnerStatus, CommandLineJobRunner}
import org.ggf.drmaa._
import java.util.{Date, Collections}

/**
 * Runs jobs using DRMAA.
 */
class DrmaaJobRunner(val session: Session, val function: CommandLineFunction) extends CommandLineJobRunner with Logging {
  /** Job Id of the currently executing job. */
  var jobId: String = _
  override def jobIdString = jobId

  // Set the display name to < 512 characters of the description
  // NOTE: Not sure if this is configuration specific?
  protected val jobNameLength = 500
  protected val jobNameFilter = """[^A-Za-z0-9_]"""
  protected def functionNativeSpec = function.jobNativeArgs.mkString(" ")

  def start() {
    session.synchronized {
      val drmaaJob: JobTemplate = session.createJobTemplate

      drmaaJob.setJobName(function.jobRunnerJobName.take(jobNameLength).replaceAll(jobNameFilter, "_"))

      // Set the current working directory
      drmaaJob.setWorkingDirectory(function.commandDirectory.getPath)

      // Set the output file for stdout
      drmaaJob.setOutputPath(":" + function.jobOutputFile.getPath)

      // If the error file is set specify the separate output for stderr
      // Otherwise join with stdout
      if (function.jobErrorFile != null) {
        drmaaJob.setErrorPath(":" + function.jobErrorFile.getPath)
      } else {
        drmaaJob.setJoinFiles(true)
      }

      if(!function.wallTime.isEmpty)
    	  drmaaJob.setHardWallclockTimeLimit(function.wallTime.get)      
      
      drmaaJob.setNativeSpecification(functionNativeSpec)

      // Instead of running the function.commandLine, run "sh <jobScript>"
      drmaaJob.setRemoteCommand("sh")
      drmaaJob.setArgs(Collections.singletonList(jobScript.toString))

      // Allow advanced users to update the request via QFunction.updateJobRun()
      updateJobRun(drmaaJob)

      updateStatus(RunnerStatus.RUNNING)

      // Start the job and store the id so it can be killed in tryStop
      try {
        Retry.attempt(() => {
          try {
            jobId = session.runJob(drmaaJob)
          } catch {
            case de: DrmaaException => throw new QException("Unable to submit job: " + de.getLocalizedMessage)
          }
        }, 1, 5, 10)
      } finally {
        // Prevent memory leaks
        session.deleteJobTemplate(drmaaJob)
      }
      logger.info("Submitted job id: " + jobId)
    }
  }

  def updateJobStatus() = {
    session.synchronized {
      var returnStatus: RunnerStatus.Value = null

      try {
        val jobStatus = session.getJobProgramStatus(jobId);
        jobStatus match {
          case Session.QUEUED_ACTIVE => returnStatus = RunnerStatus.RUNNING
          case Session.DONE =>
            val jobInfo: JobInfo = session.wait(jobId, Session.TIMEOUT_NO_WAIT)

            // Update jobInfo
            def convertDRMAATime(key: String): Date = {
              val v = jobInfo.getResourceUsage.get(key)
              if ( v != null ) new Date(v.toString.toDouble.toLong * 1000) else null;
            }
            if ( jobInfo.getResourceUsage != null ) {
              getRunInfo.startTime = convertDRMAATime("start_time")
              getRunInfo.doneTime = convertDRMAATime("end_time")
              getRunInfo.exechosts = "unknown"
            }

            if ((jobInfo.hasExited && jobInfo.getExitStatus != 0)
                || jobInfo.hasSignaled
                || jobInfo.wasAborted)
              returnStatus = RunnerStatus.FAILED
            else
              returnStatus = RunnerStatus.DONE
          case Session.FAILED => returnStatus = RunnerStatus.FAILED
          case Session.UNDETERMINED => logger.warn("Unable to determine status of job id " + jobId)
          case _ => returnStatus = RunnerStatus.RUNNING
        }
      } catch {
        // getJobProgramStatus will throw an exception once wait has run, as the
        // job will be reaped.  If the status is currently DONE or FAILED, return
        // the status.
        case de: DrmaaException =>
          if (lastStatus == RunnerStatus.DONE || lastStatus == RunnerStatus.FAILED)
            returnStatus = lastStatus
          else
            logger.warn("Unable to determine status of job id " + jobId, de)
      }

      if (returnStatus != null) {
        updateStatus(returnStatus)
        true
      } else {
        false
      }
    }
  }

  def tryStop() {
    session.synchronized {
      // Assumes that after being set the job may be
      // reassigned but will not be reset back to null
      if (jobId != null) {
        try {
          // Stop runners. SIGTERM(15) is preferred to SIGKILL(9).
          // Only way to send SIGTERM is for the Sys Admin set the terminate_method
          // resource of the designated queue to SIGTERM
          session.control(jobId, Session.TERMINATE)
        } catch {
          case e: Exception =>
            logger.error("Unable to kill job " + jobId, e)
        }
      }
    }
  }
}
