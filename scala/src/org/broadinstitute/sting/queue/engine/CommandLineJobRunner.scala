package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.IOUtils
import java.io.File

/**
 * Runs a command line function.
 */
trait CommandLineJobRunner extends JobRunner[CommandLineFunction] {
  /** A generated exec shell script. */
  protected var exec: File = _

  /** Which directory to use for the job status files. */
  protected def jobStatusDir = function.jobTempDir

  /**
   * Writes the function command line to an exec file.
   */
  protected def writeExec() {
    this.exec = IOUtils.writeTempFile(function.commandLine, ".exec", "", jobStatusDir)
  }

  /**
   * Removes all temporary files used for this LSF job.
   */
  def removeTemporaryFiles() = {
    IOUtils.tryDelete(exec)
  }

  /**
   * Builds a command line that can be run to force an automount of the directories.
   * @param function Function to look jobDirectories.
   * @return A "cd '<dir_1>' [&& cd '<dir_n>']" command.
   */
  protected def mountCommand(function: CommandLineFunction) = {
    var dirs = Set.empty[File]
    for (dir <- function.jobDirectories)
      dirs += IOUtils.dirLevel(dir, 2)
    if (dirs.size > 0)
      Some(dirs.mkString("cd '", "' && cd '", "'"))
    else
      None
  }
}
