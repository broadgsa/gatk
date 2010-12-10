package org.broadinstitute.sting.queue.engine

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util._

/**
 * Runs jobs on an LSF compute cluster.
 */
abstract class LsfJobRunner(val function: CommandLineFunction) extends DispatchJobRunner with Logging {
  protected var runStatus: RunnerStatus.Value = _

  var jobId = -1L

  /** Which directory to use for the job status files. */
  protected def jobStatusDir = function.jobTempDir

  /** A file to look for to validate that the function ran to completion. */
  protected var jobStatusPath: String = _

  /** A temporary job done file to let Queue know that the process ran successfully. */
  private lazy val jobDoneFile = new File(jobStatusPath + ".done")

  /** A temporary job done file to let Queue know that the process exited with an error. */
  private lazy val jobFailFile = new File(jobStatusPath + ".fail")

  /** A generated exec shell script. */
  protected var exec: File = _

  /** A generated pre-exec shell script. */
  protected var preExec: File = _

  /** A generated post-exec shell script. */
  protected var postExec: File = _

  // TODO: Full bsub command for debugging.
  protected def bsubCommand = "bsub " + function.commandLine

  /**
   * Updates and returns the status by looking for job status files.
   * After the job status files are detected they are cleaned up from
   * the file system and the status is cached.
   *
   * Note, these temporary job status files are currently different from the
   * .done files used to determine if a file has been created successfully.
   */
  def status = {
    try {
      if (logger.isDebugEnabled) {
        logger.debug("Done %s exists = %s".format(jobDoneFile, jobDoneFile.exists))
        logger.debug("Fail %s exists = %s".format(jobFailFile, jobFailFile.exists))
      }

      if (jobFailFile.exists) {
        removeTemporaryFiles()
        runStatus = RunnerStatus.FAILED
        logger.info("Error: " + bsubCommand)
        tailError()
      } else if (jobDoneFile.exists) {
        removeTemporaryFiles()
        runStatus = RunnerStatus.DONE
        logger.info("Done: " + bsubCommand)
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
   * Removes all temporary files used for this LSF job.
   */
  def removeTemporaryFiles() = {
    IOUtils.tryDelete(exec)
    IOUtils.tryDelete(preExec)
    IOUtils.tryDelete(postExec)
    IOUtils.tryDelete(jobDoneFile)
    IOUtils.tryDelete(jobFailFile)
  }

  /**
   * Outputs the last lines of the error logs.
   */
  protected def tailError() = {
    val errorFile = if (function.jobErrorFile != null) function.jobErrorFile else function.jobOutputFile
    if (IOUtils.waitFor(errorFile, 120)) {
      val tailLines = IOUtils.tail(errorFile, 100)
      val nl = "%n".format()
      logger.error("Last %d lines of %s:%n%s".format(tailLines.size, errorFile, tailLines.mkString(nl)))
    } else {
      logger.error("Unable to access log file: %s".format(errorFile))
    }
  }

  /**
   * Writes an exec file to cleanup any status files and
   * optionally mount any automount directories on the node.
   * @return the file path to the pre-exec.
   */
  protected def writeExec() = {
    IOUtils.writeTempFile(function.commandLine, ".exec", "", jobStatusDir)
  }

  /**
   * Writes a pre-exec file to cleanup any status files and
   * optionally mount any automount directories on the node.
   * @return the file path to the pre-exec.
   */
  protected def writePreExec() = {
    val preExec = new StringBuilder

    preExec.append("rm -f '%s/'.$LSB_JOBID.done%n".format(jobStatusDir))
    function.doneOutputs.foreach(file => preExec.append("rm -f '%s'%n".format(file)))
    preExec.append("rm -f '%s/'.$LSB_JOBID.fail%n".format(jobStatusDir))
    function.failOutputs.foreach(file => preExec.append("rm -f '%s'%n".format(file)))

    mountCommand(function).foreach(command =>
      preExec.append("%s%n".format(command)))

    IOUtils.writeTempFile(preExec.toString, ".preExec", "", jobStatusDir)
  }

  /**
   * Writes a post-exec file to create the status files.
   * @return the file path to the post-exec.
   */
  protected def writePostExec() = {
    val postExec = new StringBuilder

    val touchDone = function.doneOutputs.map("touch '%s'%n".format(_)).mkString
    val touchFail = function.failOutputs.map("touch '%s'%n".format(_)).mkString

    postExec.append("""|
  |if [ "${LSB_JOBPEND:-unset}" != "unset" ]; then
  |  exit 0
  |fi
  |
  |JOB_STAT_ROOT='%s/'.$LSB_JOBID
  |if [ "$LSB_JOBEXIT_STAT" == "0" ]; then
  |%stouch "$JOB_STAT_ROOT".done
  |else
  |%stouch "$JOB_STAT_ROOT".fail
  |fi
  |""".stripMargin.format(jobStatusDir, touchDone, touchFail))

    IOUtils.writeTempFile(postExec.toString, ".postExec", "", jobStatusDir)
  }
}
