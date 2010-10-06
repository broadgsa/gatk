package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.Input
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.apache.commons.io.FileUtils

/**
 * Removes the temporary directories for scatter / gather.
 * The script can be changed by setting rmdirScript.
 * By default uses rm -rf.
 * The format of the call is <mkdirScript> <dir_1> [.. <dir_n>]
 */
class CleanupTempDirsFunction extends InProcessFunction {
  @Input(doc="Original outputs of the gather functions")
  var originalOutputs: Set[File] = Set.empty[File]

  @Input(doc="Temporary directories to be deleted")
  var tempDirectories: List[File] = Nil

  def run() = tempDirectories.foreach(FileUtils.deleteDirectory(_))
}
