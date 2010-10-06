package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.QException

/**
 * Runs a job on the command line by invoking "sh -c <command>"
 */
class ShellJob extends CommandLineJob with Logging {
  /**
   * Runs the command and waits for the output.
   */
  def run() = {
    assert(command != null, "Command was not set on job")

    val (redirectError, errorFile) = if (this.errorFile == null) (true, null) else (false, this.errorFile)
    val bufferSize = if (redirectError || logger.isDebugEnabled) FIVE_MB else 0
    val stdinSettings = new ProcessController.InputStreamSettings(null, this.inputFile)
    val stdoutSettings = new ProcessController.OutputStreamSettings(bufferSize, this.outputFile, true)
    val stderrSettings = new ProcessController.OutputStreamSettings(FIVE_MB, errorFile, true)
    val commandLine = Array("sh", "-c", command)
    val processSettings = new ProcessController.ProcessSettings(
      commandLine, null, this.workingDir, stdinSettings, stdoutSettings, stderrSettings, redirectError)

    val output = processController.exec(processSettings)

    if (output.exitValue != 0) {
      val streamOutput = if (redirectError) output.stdout else output.stderr
      throw new JobExitException("Failed to run job.", commandLine, output.exitValue, content(streamOutput))
    }
  }
}
