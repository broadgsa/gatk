package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.commandline.{Output, Input}

class CreateTempDirsFunction extends CommandLineFunction {
  @Input(doc="Original inputs to the scattered function")
  var originalInputs: Set[Any] = Set.empty[Any]

  @Output(doc="Temporary directories to create")
  var tempDirectories: List[File] = Nil

  @Input(doc="Sleep seconds", required=false)
  var mkdirSleepSeconds: Option[Int] = None

  // TODO: After port of LSF submitter use -cwd <dir> instead of trying to run from the directory
  // For now, create the directory so that BroadCore can run bsub from it -kshakir July 27, 2010 on chartl's computer

  override def freeze = {
    super.freeze
    tempDirectories.foreach(_.mkdirs)
  }

  def commandLine = "mkdir -pv%s%s".format(repeat(" '", tempDirectories, "'"), optional(" && sleep ", mkdirSleepSeconds))
}
