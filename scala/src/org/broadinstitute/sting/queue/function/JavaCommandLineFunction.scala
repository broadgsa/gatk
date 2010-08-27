package org.broadinstitute.sting.queue.function

/**
 * Defines a command line function that runs java code.
 */
trait JavaCommandLineFunction extends CommandLineFunction {
  /**
   * Returns the java executable to run.
   */
  def javaExecutable: String

  def javaOpts = "%s -Djava.io.tmpdir=%s"
    .format(optional(" -Xmx", memoryLimit, "g"), jobTempDir)

  def commandLine = "java%s %s"
    .format(javaOpts, javaExecutable)
}
