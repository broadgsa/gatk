/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.function

import org.broadinstitute.gatk.queue.util._
import org.broadinstitute.gatk.utils.commandline.Argument

/**
 * A command line that will be run in a pipeline.
 */
trait CommandLineFunction extends QFunction with Logging {
  def commandLine: String

  /** Setting the wall time request for DRMAA / run limit for LSF */
  var wallTime: Option[Long] = None
  
  /** Upper memory limit */
  @Argument(doc="Memory limit", required=false)
  var memoryLimit: Option[Double] = None

  /** Resident memory limit */
  @Argument(doc="Resident memory limit", required=false)
  var residentLimit: Option[Double] = None

  /** Resident memory request */
  @Argument(doc="Resident memory request", required=false)
  var residentRequest: Option[Double] = None

  /** the number of SMP cores this job wants */
  var nCoresRequest: Option[Int] = None

  /** Job project to run the command */
  var jobProject: String = _

  /** Job queue to run the command */
  var jobQueue: String = _

  /** Native arguments to pass to the job runner */
  var jobNativeArgs: Seq[String] = Nil

  /** Native arguments to pass to the job runner */
  var jobResourceRequests: Seq[String] = Nil

  /** Environment names to pass to the job runner */
  var jobEnvironmentNames: Seq[String] = Nil

  override def copySettingsTo(function: QFunction) {
    super.copySettingsTo(function)
    function match {
      case commandLineFunction: CommandLineFunction =>
        if(commandLineFunction.wallTime.isEmpty)
          commandLineFunction.wallTime = this.wallTime
        
        if (commandLineFunction.memoryLimit.isEmpty)
          commandLineFunction.memoryLimit = this.memoryLimit

        if (commandLineFunction.residentLimit.isEmpty)
          commandLineFunction.residentLimit = this.residentLimit

        if (commandLineFunction.residentRequest.isEmpty)
          commandLineFunction.residentRequest = this.residentRequest

        if (commandLineFunction.nCoresRequest.isEmpty)
          commandLineFunction.nCoresRequest = this.nCoresRequest

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
   
    if(wallTime.isEmpty)
      wallTime = qSettings.jobWalltime
    
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

    // the default value is 1 core
    if (nCoresRequest.isEmpty)
      nCoresRequest = Some(1)

    if (residentRequest.isEmpty)
      residentRequest = memoryLimit

    if (residentLimit.isEmpty || residentLimit == residentRequest)
      residentLimit = residentRequest.map(residentLimitBuffer)

    super.freezeFieldValues()
  }

  /**
   * @return A function that decides how much memory cushion to add to the residentRequest to create the residentLimit
   */
  def residentLimitBuffer: (Double => Double) = (1.2 * _)

  /**
   * Safely construct a full required command-line argument with consistent quoting, whitespace separation, etc.
   *
   * @param prefix Prefix to insert before the argument value (eg., "-f")
   * @param param The argument value itself
   * @param suffix Suffix to append after the argument value
   * @param spaceSeparated If true, insert a space between the prefix, param, and suffix
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @param format Format String used to convert param to a String
   * @return The combined and formatted argument, surrounded by whitespace
   */
  protected def required( prefix: String, param: Any, suffix: String = "", spaceSeparated: Boolean = true,
                          escape: Boolean = true, format: String = "%s" ): String = {
    " %s ".format(formatArgument(prefix, param, suffix, spaceSeparated, escape, format))
  }

  /**
   * Safely construct a one-token required command-line argument with quoting
   *
   * @param param The command-line argument value
   * @return The argument value quoted and surrounded by whitespace
   */
  protected def required( param: Any ): String = {
    required("", param)
  }

  /**
   * Safely construct a one-token required command-line argument, and specify whether you want quoting
   *
   * @param param The command-line argument value
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @return The argument value, quoted if quoting was requested, and surrounded by whitespace
   */
  protected def required( param: Any, escape: Boolean ): String = {
    required("", param, escape=escape)
  }

  /**
   * Safely construct a full optional command-line argument with consistent quoting, whitespace separation, etc.
   * If the argument has no value, returns an empty String.
   *
   * @param prefix Prefix to insert before the argument value (eg., "-f")
   * @param param The argument value itself (if null/empty, the method returns empty String)
   * @param suffix Suffix to append after the argument value
   * @param spaceSeparated If true, insert a space between the prefix, param, and suffix
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @param format Format String used to convert param to a String
   * @return The combined and formatted argument, surrounded by whitespace, or an empty String
   *         if the argument has no value
   */
  protected def optional( prefix: String, param: Any, suffix: String = "", spaceSeparated: Boolean = true,
                          escape: Boolean = true, format: String = "%s" ): String = {
    if ( hasValue(param) ) " %s ".format(formatArgument(prefix, param, suffix, spaceSeparated, escape, format)) else ""
  }

  /**
   * Safely construct a one-token optional command-line argument with quoting.
   * If the argument has no value, returns an empty String.
   *
   * @param param The command-line argument value
   * @return The argument value quoted and surrounded by whitespace, or an empty String
   *         if the argument has no value
   */
  protected def optional( param: Any ): String = {
    optional("", param)
  }

  /**
   * Safely construct a one-token conditional command-line argument. If the provided condition
   * is false, an empty String is returned.
   *
   * @param condition The condition to check
   * @param param The command-line argument value
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @param format Format String used to convert param to a String
   * @return The command-line argument value, quoted if quoting was requested and surrounded
   *         by whitespace, or an empty String if the argument has no value.
   */
  protected def conditional( condition: Boolean, param: Any, escape: Boolean = true, format: String = "%s" ): String = {
    if ( condition ) {
      " %s ".format(formatArgument("", param, "", spaceSeparated = false, escape = escape, paramFormat = format))
    }
    else {
      ""
    }
  }

  /**
   * Safely construct a series of full command-line arguments with consistent quoting, whitespace separation, etc.
   *
   * Each argument value is preceded by a prefix/suffix if they are set. A function can be provided to vary
   * each prefix for each argument value (eg., -f:tag1 file1 -f:tag2 file2) -- the default is to use
   * the same prefix for all arguments.
   *
   * @param prefix Prefix to insert before each argument value (eg., "-f")
   * @param params The collection of argument values
   * @param suffix Suffix to append after each argument value
   * @param separator Specifies how to separate the various arguments from each other
   *                  (eg., what should go between '-f' 'file1' and '-f' 'file2'?)
   *                  Default is one space character.
   * @param spaceSeparated If true, insert a space between each individual prefix, param, and suffix
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @param format Format String used to convert each individual param within params to a String
   * @param formatPrefix Function mapping (prefix, argumentValue) pairs to prefixes. Can be used to
   *                     vary each prefix depending on the argument value (useful for tags, etc.).
   *                     Default is to use the same prefix for all argument values.
   * @return The series of command-line arguments, quoted and whitespace-delimited as requested,
   *         or an empty String if params was null/Nil/None.
   */
  protected def repeat(prefix: String, params: Traversable[_], suffix: String = "", separator: String = " ",
                       spaceSeparated: Boolean = true, escape: Boolean = true, format: String = "%s",
                       formatPrefix: (String, Any) => String = (prefix, value) => prefix): String = {
    if (CollectionUtils.isNullOrEmpty(params))
      ""
    else
      " %s ".format(params.filter(param => hasValue(param)).map(param => formatArgument(formatPrefix(prefix, param), param, suffix, spaceSeparated, escape, format)).mkString(separator))
  }

  /**
   * Safely construct a series of one-token command-line arguments with quoting and space separation.
   *
   * @param params The collection of argument values
   * @return The argument values quoted and space-delimited, or an empty String if params was null/Nil/None
   */
  protected def repeat( params: Traversable[_] ): String = {
    repeat("", params)
  }

  /**
   * Given an (optional) prefix, an argument value, and an (optional) suffix, formats a command-line
   * argument with the specified level of quoting and space-separation.
   *
   * Helper method for required(), optional(), conditional(), and repeat() -- do not use this
   * method directly!
   *
   * @param prefix Prefix to insert before the argument value (eg., "-f"). Ignored if empty/null.
   * @param param The argument value itself. If this is Some(x), it is unwrapped to x before processing.
   * @param suffix Suffix to append after the argument value. Ignored if empty/null.
   * @param spaceSeparated If true, insert a space between the prefix, param, and suffix
   * @param escape If true, quote the generated argument to avoid interpretation by the shell
   * @param paramFormat Format string used to convert param to a String
   * @return The combined and formatted argument, NOT surrounded by any whitespace.
   *         Returns an empty String if param was null/empty.
   */
  protected def formatArgument( prefix: String, param: Any, suffix: String, spaceSeparated: Boolean, escape: Boolean,
                                paramFormat: String ): String = {
    if (CollectionUtils.isNullOrEmpty(param)) {
      return ""
    }

    // Trim leading and trailing whitespace off our three tokens, and unwrap Some(x) to x for the param
    val trimmedValues : Seq[String] = Seq((if ( prefix != null ) prefix.trim else ""),
                                            (param match {
                                              case Some(x) => paramFormat.format(x).trim
                                              case x => paramFormat.format(x).trim
                                            }),
                                            (if ( suffix != null ) suffix.trim else ""))
    var joinedArgument : String = null

    // If the user requested space-separation, join the tokens with a space, and escape individual
    // NON-EMPTY tokens if escaping was requested (eg., ("-f", "foo", "") -> "'-f' 'foo'")
    if ( spaceSeparated ) {
      joinedArgument = trimmedValues.map(x => if ( x.length > 0 && escape ) ShellUtils.escapeShellArgument(x) else x).mkString(" ").trim()
    }

    // Otherwise join the tokens without any intervening whitespace, and if quoting was requested
    // quote the entire concatenated value (eg., ("-Xmx", "4", "G") -> "'-Xmx4G'")
    else  {
      joinedArgument = if ( escape ) ShellUtils.escapeShellArgument(trimmedValues.mkString("")) else trimmedValues.mkString("")
    }

    // If the user requested escaping and we ended up with an empty String after joining, quote the empty
    // String to preserve the command line token. Otherwise just return the joined argument
    if ( joinedArgument.length == 0 && escape ) ShellUtils.escapeShellArgument(joinedArgument) else joinedArgument
  }
}
