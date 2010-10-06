package org.broadinstitute.sting.queue.util

/**
 * bkills a list of lsf jobs.
 */
class LsfKillJob(jobs: List[LsfJob]) extends CommandLineJob with Logging {
  command = "bkill " + jobs.map(_.bsubJobId).mkString(" ")

  def run() = {
    // capture the output for debugging
    val stdinSettings = new ProcessController.InputStreamSettings(null, null)
    val stdoutSettings = new ProcessController.OutputStreamSettings(FIVE_MB, null, false)
    val stderrSettings = new ProcessController.OutputStreamSettings(FIVE_MB, null, false)

    val bkillCommand = (List("bkill") ++ jobs.map(_.bsubJobId)).toArray

    // launch the bsub job from the current directory
    val processSettings = new ProcessController.ProcessSettings(
      bkillCommand, null, null, stdinSettings, stdoutSettings, stderrSettings, false)
    val bkillOutput = processController.exec(processSettings)

    if (bkillOutput.exitValue != 0) {
      throw new JobExitException("Failed to kill LSF jobs.", bkillCommand, bkillOutput.exitValue, content(bkillOutput.stderr))
    }
  }
}
