package org.broadinstitute.sting.queue.function

/**
 * Defines a command line function that runs java code.
 */
trait JavaCommandLineFunction extends CommandLineFunction {
  /**
   * Returns the java executable to run.
   */
  def javaExecutable: String

  /**
   * Memory limit for the java executable, or if None will use the default memoryLimit.
   */
  var javaMemoryLimit: Option[Int] = None

  override def freezeFieldValues = {
    super.freezeFieldValues

    if (javaMemoryLimit.isEmpty && memoryLimit.isDefined)
      javaMemoryLimit = memoryLimit
  }

  def javaOpts = "%s -Djava.io.tmpdir=%s"
    .format(optional(" -Xmx", javaMemoryLimit, "g"), jobTempDir)

  def commandLine = "java%s %s"
    .format(javaOpts, javaExecutable)
}
