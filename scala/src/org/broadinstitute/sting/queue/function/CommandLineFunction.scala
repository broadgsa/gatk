package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.commandline._
import java.io.File
import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.scattergather.{SimpleTextGatherFunction, Gather}

/**
 * A command line that will be run in a pipeline.
 */
trait CommandLineFunction extends QFunction with Logging {
  def commandLine: String

  /** Upper memory limit */
  var memoryLimit: Option[Int] = None

  /** Whether a job is restartable */
  var jobRestartable = true

  /** Prefix for automatic job name creation */
  var jobNamePrefix: String = _

  /** The name name of the job */
  var jobName: String = _

  /** Job project to run the command */
  var jobProject: String = _

  /** Job queue to run the command */
  var jobQueue: String = _

  /** Temporary directory to write any files */
  var jobTempDir: File = IOUtils.javaTempDir

  /** File to redirect any output.  Defaults to <jobName>.out */
  @Output(doc="File to redirect any output", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobOutputFile: File = _

  /** File to redirect any errors.  Defaults to <jobName>.out */
  @Output(doc="File to redirect any errors", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobErrorFile: File = _

  /**
   * Returns set of directories required to run the command.
   * @return Set of directories required to run the command.
   */
  def jobDirectories = {
    var dirs = Set.empty[File]
    dirs += commandDirectory
    if (jobTempDir != null)
      dirs += jobTempDir
    dirs ++= inputs.map(_.getParentFile)
    dirs ++= outputs.map(_.getParentFile)
    dirs
  }

  override protected def useStatusOutput(file: File) =
    file != jobOutputFile && file != jobErrorFile

  override def description = commandLine

  /**
   * The function description in .dot files
   */
  override def dotString = jobName + " => " + commandLine

  /**
   * Sets all field values.
   */
  override def freezeFieldValues = {
    if (jobNamePrefix == null)
      jobNamePrefix = qSettings.jobNamePrefix

    if (jobQueue == null)
      jobQueue = qSettings.jobQueue

    if (jobProject == null)
      jobProject = qSettings.jobProject

    if (memoryLimit.isEmpty && qSettings.memoryLimit.isDefined)
      memoryLimit = qSettings.memoryLimit

    if (jobName == null)
      jobName = CommandLineFunction.nextJobName(jobNamePrefix)

    if (jobOutputFile == null)
      jobOutputFile = new File(jobName + ".out")

    super.freezeFieldValues
  }

  /**
   * Repeats parameters with a prefix/suffix if they are set otherwise returns "".
   * Skips null, Nil, None.  Unwraps Some(x) to x.  Everything else is called with x.toString.
   * @param prefix Command line prefix per parameter.
   * @param params Traversable parameters.
   * @param suffix Optional suffix per parameter.
   * @param separator Optional separator per parameter.
   * @param format Format function if the value has a value
   * @return The generated string
   */
  protected def repeat(prefix: String, params: Traversable[_], suffix: String = "", separator: String = "",
                       format: (String, Any, String) => String = formatValue("%s")) =
    params.filter(param => hasValue(param)).map(param => format(prefix, param, suffix)).mkString(separator)

  /**
   * Returns parameter with a prefix/suffix if it is set otherwise returns "".
   * Does not output null, Nil, None.  Unwraps Some(x) to x.  Everything else is called with x.toString.
   * @param prefix Command line prefix per parameter.
   * @param param Parameter to check for a value.
   * @param suffix Optional suffix per parameter.
   * @param format Format function if the value has a value
   * @return The generated string
   */
  protected def optional(prefix: String, param: Any, suffix: String = "",
                         format: (String, Any, String) => String = formatValue("%s")) =
    if (hasValue(param)) format(prefix, param, suffix) else ""

  /**
   * Returns "" if the value is null or an empty collection, otherwise return the value.toString.
   * @param format Format string if the value has a value
   * @param prefix Command line prefix per parameter.
   * @param param Parameter to check for a value.
   * @param suffix Optional suffix per parameter.
   * @return "" if the value is null, or "" if the collection is empty, otherwise the value.toString.
   */
  protected def formatValue(format: String)(prefix: String, param: Any, suffix: String): String =
    if (CollectionUtils.isNullOrEmpty(param))
      ""
    else
      prefix + (param match {
        case Some(x) => format.format(x)
        case x => format.format(x)
      }) + suffix

}

/**
 * A command line that will be run in a pipeline.
 */
object CommandLineFunction {
  /** Job index counter for this run of Queue. */
  private var jobIndex = 0

  /**
   * Returns the next job name using the prefix.
   * @param prefix Prefix of the job name.
   * @return the next job name.
   */
  private def nextJobName(prefix: String) = {
    jobIndex += 1
    prefix + "-" + jobIndex
  }
}
