package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.FileProvider

/**
 * Used to provide -I input_file arguments to the GATK.
 */
class NamedFile(var file: File, var name: String = null) extends FileProvider {
  require(file != null, "NamedFile file cannot be null")
}

/**
 * Used to provide -I input_file arguments to the GATK.
 */
object NamedFile {
  /**
   * Formats the rod binding on the command line.
   * Used for optional and repeat.
   * @param cmdLineParam command line parameter, ex: -I
   * @param prefix unused
   * @param value NamedFile to add.
   * @param suffix unused
   * @return The command line addition.
   */
  def formatCommandLine(cmdLineParam: String)(prefix: String, value: Any, suffix: String) = {
    value match {
      case namedFile: NamedFile =>
        if (namedFile.name != null)
          " %s:%s %s".format(cmdLineParam, namedFile.name, namedFile.file)
        else
          " %s %s".format(cmdLineParam, namedFile.file)
    }
  }
}
