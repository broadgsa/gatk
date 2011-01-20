package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.queue.util.{Logging, IOUtils}

/**
 * Runs a command line function.
 */
trait CommandLineJobRunner extends JobRunner[CommandLineFunction] with Logging {
  /** A generated exec shell script. */
  protected var exec: File = _

  /** Which directory to use for the job status files. */
  protected def jobStatusDir = function.jobTempDir

  /** The last time the status returned unknown. */
  protected var firstUnknownTime: Option[Long] = None

  /** Amount of time the status can return unknown before giving up. */
  protected val unknownStatusMaxSeconds = 5 * 60

  /** Number of seconds for a non-normal exit status before we give up on expecting LSF to retry the function. */
  protected val retryExpiredSeconds = 5 * 60

  /**
   * Writes the function command line to an exec file.
   */
  protected def writeExec() {
    var exec = new StringBuilder
    
    var dirs = function.mountDirectories
    for (dir <- function.jobDirectories)
      dirs += IOUtils.dirLevel(dir, 2)
    if (dirs.size > 0) {
      // prepend "cd '<dir_1>' [&& cd '<dir_n>']" to automount the directories.
      exec.append(dirs.mkString("cd '", "' && cd '", "'"))
      exec.append(" && cd '%s' && \\%n".format(function.commandDirectory))
    }
    exec.append(function.commandLine)

    this.exec = IOUtils.writeTempFile(exec.toString, ".exec", "", jobStatusDir)
  }

  /**
   * Removes all temporary files used for this LSF job.
   */
  def removeTemporaryFiles() = {
    IOUtils.tryDelete(exec)
  }

  /**
   * Outputs the last lines of the error logs.
   */
  protected def tailError() = {
    val errorFile = functionErrorFile
    if (IOUtils.waitFor(errorFile, 120)) {
      val tailLines = IOUtils.tail(errorFile, 100)
      val nl = "%n".format()
      logger.error("Last %d lines of %s:%n%s".format(tailLines.size, errorFile, tailLines.mkString(nl)))
    } else {
      logger.error("Unable to access log file: %s".format(errorFile))
    }
  }
}
