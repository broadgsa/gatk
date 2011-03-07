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

  /**
   * Writes the function command line to an exec file.
   */
  protected def writeExec() {
    var exec = new StringBuilder
    
    var dirs = Set.empty[File]
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

  override def removeTemporaryFiles() {
    super.removeTemporaryFiles()
    IOUtils.tryDelete(exec)
  }
}
