package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.commandline.{Argument, Input}

/**
 * Removes the temporary directories for scatter / gather.
 * The script can be changed by setting rmdirScript.
 * By default uses rm -rf.
 * The format of the call is <mkdirScript> <dir_1> [.. <dir_n>]
 */
class CleanupTempDirsFunction extends CommandLineFunction {
  @Input(doc="Original outputs of the gather functions")
  var originalOutputs: Set[File] = Set.empty[File]

  @Input(doc="Temporary directories to be deleted")
  var tempDirectories: List[File] = Nil

  @Argument(doc="rmdir script or command")
  var rmdirScript = "rm -rf"

  def commandLine = "%s%s".format(rmdirScript, repeat(" '", tempDirectories, "'"))
}
