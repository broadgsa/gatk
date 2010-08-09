package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.commandline.Argument
import java.io.File

/**
 * Defines a command line function that runs from a jar file.
 */
trait JarCommandLineFunction extends CommandLineFunction {
  @Argument(doc="jar")
  var jarFile: File = _

  def commandLine = "java%s -Djava.io.tmpdir=%s -jar %s"
    .format(optional(" -Xmx", memoryLimit, "g"), jobTempDir, jarFile)
}
