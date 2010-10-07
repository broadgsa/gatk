package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{Output, Input}
import org.broadinstitute.sting.queue.function.InProcessFunction

/**
 * Creates the temporary directories for scatter / gather.
 * The script can be changed by setting mkdirScript.
 * By default uses mkdir -pv
 * The format of the call is <rmdirScript> <dir_1> [.. <dir_n>]
 */
class CreateTempDirsFunction extends InProcessFunction {
  @Input(doc="Original inputs to the scattered function")
  var originalInputs: Set[File] = Set.empty[File]

  @Output(doc="Temporary directories to create")
  var tempDirectories: List[File] = Nil

  override def useStatusOutput(file: File) = false

  def run() = tempDirectories.foreach(_.mkdirs)
}
