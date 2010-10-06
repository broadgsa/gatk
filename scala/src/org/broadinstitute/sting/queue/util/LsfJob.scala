package org.broadinstitute.sting.queue.util

import java.util.regex.Pattern
import collection.JavaConversions._

/**
 * A job submitted to LSF. This class is designed to work somewhat like
 * java.lang.Process, but has some extensions.
 *
 * @author A subset of the original BroadCore ported to scala by Khalid Shakir
 */
class LsfJob extends CommandLineJob with Logging {
  var name: String = _
  var project: String = _
  var queue: String = _
  var preExecCommand: String = _
  var postExecCommand: String = _
  var waitForCompletion = false
  var extraBsubArgs: List[String] = Nil
  var bsubJobId: String = _

  /**
   * Starts the job. Command must exist. The job will be submitted to LSF.
   */
  def run() = {
    assert(bsubJobId == null, "LSF job was already submitted")
    assert(command != null, "Command was not set on LSF job")
    assert(outputFile != null, "Output file must be set on LSF job")

    // capture the output for debugging
    val stdinSettings = new ProcessController.InputStreamSettings(null, null)
    val stdoutSettings = new ProcessController.OutputStreamSettings(FIVE_MB, null, false)
    val stderrSettings = new ProcessController.OutputStreamSettings(FIVE_MB, null, false)

    // launch the bsub job from the current directory
    val processSettings = new ProcessController.ProcessSettings(
      bsubCommand, environmentVariables, null, stdinSettings, stdoutSettings, stderrSettings, false)
    val bsubOutput = processController.exec(processSettings)

    if (bsubOutput.exitValue != 0) {
      throw new JobExitException("Failed to submit LSF job.", bsubCommand,
        bsubOutput.exitValue, content(bsubOutput.stderr))
    }

    // get the LSF job ID
    val matcher = LsfJob.JOB_ID.matcher(bsubOutput.stdout.content)
    matcher.find()
    bsubJobId = matcher.group
  }

  /**
   * Generates the bsub command line for this LsfJob.
   * @return command line as a Array[String]
   */
  def bsubCommand = {
    var args = List.empty[String]
    args :+= "bsub"

    if (name != null) {
      args :+= "-J"
      args :+= name
    }

    if (inputFile != null) {
      args :+= "-i"
      args :+= inputFile.getAbsolutePath
    }

    args :+= "-o"
    args :+= outputFile.getAbsolutePath

    if (errorFile != null) {
      args :+= "-e"
      args :+= errorFile.getAbsolutePath
    }

    if (queue != null) {
      args :+= "-q"
      args :+= queue
    }

    if (project != null) {
      args :+= "-P"
      args :+= project
    }

    if (preExecCommand != null) {
      args :+= "-E"
      args :+= preExecCommand
    }

    if (postExecCommand != null) {
      args :+= "-Ep"
      args :+= postExecCommand
    }

    if (workingDir != null) {
      args :+= "-cwd"
      args :+= workingDir.getPath
    }

    if (waitForCompletion) {
      args :+= "-K"
    }

    args ++= extraBsubArgs

    args :+= command

    args.toArray
  }

  /**
   * Get the list of environment variables and pass into the exec job.  We strip
   * out LD_ASSUME_KERNEL because it behaves badly when running bsub jobs across
   * different versions of the linux OS.
   *
   * @return array of environment vars in 'name=value' format.
   */
  private def environmentVariables =
    System.getenv()
            .filterNot{case (name, value) => name.equalsIgnoreCase("LD_ASSUME_KERNEL") || value == null}
            .toMap
}

/**
 * A job submitted to LSF. This class is designed to work somewhat like
 * java.lang.Process, but has some extensions.
 *
 * @author A subset of the original BroadCore ported to scala by Khalid Shakir
 */
object LsfJob {
  /** Used to search the stdout for the job id. */
  private val JOB_ID = Pattern.compile("\\d+")
}
