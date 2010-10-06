package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.util.IOUtils
import java.io.File

/**
 * Dispatches jobs to a compute cluster.
 */
trait DispatchJobRunner extends JobRunner {
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
