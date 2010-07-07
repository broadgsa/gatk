package org.broadinstitute.sting.queue.function

import java.io.File
import java.lang.management.ManagementFactory
import org.broadinstitute.sting.queue.function.scattergather.{Gather, SimpleTextGatherFunction}
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.commandline.{ClassType, Output, Input}

trait DispatchFunction extends InputOutputFunction {
  def commandLine: String

  @Input(doc="Upper memory limit", required=false)
  @ClassType(classOf[Int])
  var memoryLimit: Option[Int] = None

  /**
   * The directory where the command should run.
   */
  @Input(doc="Directory to write any files", required=false)
  var commandDirectory: File = IOUtils.CURRENT_DIR

  @Input(doc="Prefix for automatic job name creation", required=false)
  var jobNamePrefix: String = _

  @Input(doc="Job name to run on the farm", required=false)
  var jobName: String = _

  @Output(doc="File to redirect any output", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobOutputFile: File = _

  @Output(doc="File to redirect any errors", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobErrorFile: File = _

  @Input(doc="Job project to run the command", required=false)
  var jobProject = "Queue"

  @Input(doc="Job queue to run the command", required=false)
  var jobQueue = "broad"

  override def freeze = {
    if (jobNamePrefix == null)
      jobNamePrefix = DispatchFunction.processNamePrefix

    if (jobName == null)
      jobName = DispatchFunction.nextJobName(jobNamePrefix)

    if (jobOutputFile == null)
      jobOutputFile = new File(jobName + ".out")

    if (jobErrorFile == null)
      jobErrorFile = new File(jobName + ".err")

    commandDirectory = IOUtils.absolute(IOUtils.CURRENT_DIR, commandDirectory)

    super.freeze
  }

  /**
   * Override the canon function to change any relative path to an absolute path.
   */
  override protected def canon(value: Any) = {
    value match {
      case file: File => IOUtils.absolute(commandDirectory, file)
      case x => super.canon(x)
    }
  }

  def absolute(file: File) = IOUtils.absolute(commandDirectory, file)
  def temp(subDir: String) = IOUtils.sub(commandDirectory, jobName + "-" + subDir)

  override def toString = commandLine
}

object DispatchFunction {
  private val processNamePrefix = "Q-" + {
    var prefix = ManagementFactory.getRuntimeMXBean.getName
    val index = prefix.indexOf(".")
    if (index >= 0)
      prefix = prefix.substring(0, index)
    prefix
  }

  private var jobIndex = 0

  private def nextJobName(prefix: String) = {
    jobIndex += 1
    prefix + "-" + jobIndex
  }
}
