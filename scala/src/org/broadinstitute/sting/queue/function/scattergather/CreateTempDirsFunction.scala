package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.commandline.{Argument, Output, Input}

/**
 * Creates the temporary directories for scatter / gather.
 * The script can be changed by setting mkdirScript.
 * By default uses mkdir -pv
 * The format of the call is <rmdirScript> <dir_1> [.. <dir_n>]
 */
class CreateTempDirsFunction extends CommandLineFunction {
  @Input(doc="Original inputs to the scattered function")
  var originalInputs: Set[File] = Set.empty[File]

  @Output(doc="Temporary directories to create")
  var tempDirectories: List[File] = Nil

  @Argument(doc="mkdir script or command")
  var mkdirScript = "mkdir -pv"

  override def upToDate = tempDirectories.forall(_.exists)

  def commandLine = "%s%s".format(mkdirScript, repeat(" '", tempDirectories, "'"))

  /**
   * This function is creating the directories, so returns just this command directory.
   */
  override def jobDirectories = Set(commandDirectory)
}
