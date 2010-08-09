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
    val bufferSize = if (logger.isDebugEnabled) FIVE_MB else 0
    val stdinSettings = new ProcessController.InputStreamSettings(null, this.inputFile)
    val stdoutSettings = new ProcessController.OutputStreamSettings(bufferSize, this.outputFile, true)
    val stderrSettings = new ProcessController.OutputStreamSettings(FIVE_MB, errorFile, true)
    val processSettings = new ProcessController.ProcessSettings(
      Array("sh", "-c", command), null, this.workingDir, stdinSettings, stdoutSettings, stderrSettings, redirectError)

    val output = processController.exec(processSettings)

    if (logger.isDebugEnabled) {
      logger.debug("output: " + content(output.stdout))
      logger.debug("error: " + content(output.stderr))
      logger.debug("Command exited with result: " + output.exitValue)
    }

    if (output.exitValue != 0) {
      logger.error("Failed to run job, got exit code %s.  Standard error contained: %n%s"
              .format(output.exitValue, content(output.stderr)))
      throw new QException("Failed to run job, got exit code %s.".format(output.exitValue))
    }
  }
}
