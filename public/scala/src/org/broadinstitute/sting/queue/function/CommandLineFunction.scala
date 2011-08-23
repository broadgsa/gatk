package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.util._

/**
 * A command line that will be run in a pipeline.
 */
trait CommandLineFunction extends QFunction with Logging {
  def commandLine: String

  /** Upper memory limit */
  var memoryLimit: Option[Double] = None

  /** Resident memory limit */
  var residentLimit: Option[Double] = None

  /** Resident memory request */
  var residentRequest: Option[Double] = None

  /** Job project to run the command */
  var jobProject: String = _

  /** Job queue to run the command */
  var jobQueue: String = _

  /** Native arguments to pass to the job runner */
  var jobNativeArgs: List[String] = Nil

  /** Native arguments to pass to the job runner */
  var jobResourceRequests: List[String] = Nil

  /** Environment names to pass to the job runner */
  var jobEnvironmentNames: List[String] = Nil

  override def copySettingsTo(function: QFunction) {
    super.copySettingsTo(function)
    function match {
      case commandLineFunction: CommandLineFunction =>
        if (commandLineFunction.memoryLimit.isEmpty)
          commandLineFunction.memoryLimit = this.memoryLimit

        if (commandLineFunction.residentLimit.isEmpty)
          commandLineFunction.residentLimit = this.residentLimit

        if (commandLineFunction.residentRequest.isEmpty)
          commandLineFunction.residentRequest = this.residentRequest

        if (commandLineFunction.jobProject == null)
          commandLineFunction.jobProject = this.jobProject

        if (commandLineFunction.jobQueue == null)
          commandLineFunction.jobQueue = this.jobQueue

        if (commandLineFunction.jobNativeArgs.isEmpty)
          commandLineFunction.jobNativeArgs = this.jobNativeArgs

        if (commandLineFunction.jobResourceRequests.isEmpty)
          commandLineFunction.jobResourceRequests = this.jobResourceRequests

        if (commandLineFunction.jobEnvironmentNames.isEmpty)
          commandLineFunction.jobEnvironmentNames = this.jobEnvironmentNames

      case _ => /* ignore */
    }
  }

  /**
   * Returns set of directories required to run the command.
   * @return Set of directories required to run the command.
   */
  def jobDirectories = outputDirectories ++ inputs.map(_.getParentFile)

  override def description = commandLine

  /**
   * Sets all field values.
   */
  override def freezeFieldValues() {
    if (jobQueue == null)
      jobQueue = qSettings.jobQueue

    if (jobProject == null)
      jobProject = qSettings.jobProject

    if (jobNativeArgs.isEmpty)
      jobNativeArgs = qSettings.jobNativeArgs

    if (jobResourceRequests.isEmpty)
      jobResourceRequests = qSettings.jobResourceRequests

    if (jobEnvironmentNames.isEmpty)
      jobEnvironmentNames = qSettings.jobEnvironmentNames

    if (memoryLimit.isEmpty)
      memoryLimit = qSettings.memoryLimit

    if (residentLimit.isEmpty)
      residentLimit = qSettings.residentLimit

    if (residentRequest.isEmpty)
      residentRequest = qSettings.residentRequest

    if (residentRequest.isEmpty)
      residentRequest = memoryLimit

    if (residentLimit.isEmpty)
      residentLimit = residentRequest.map( _ * 1.2 )

    super.freezeFieldValues()
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
    if (params == null)
      ""
    else
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

  /**
   * Returns the parameter if the condition is true. Useful for long string of parameters
   * @param condition the condition to validate
   * @param param the string to be returned in case condition is true
   * @return param if condition is true, "" otherwise
   */
  protected def conditionalParameter(condition: Boolean, param: String): String =
    if (condition == true)
      param
    else
      ""
}
