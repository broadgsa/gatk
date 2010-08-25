package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.FileProvider

/**
 * Used to provide -B rodBind arguments to the GATK.
 */
case class RodBind(var trackName: String, var trackType: String, var file: File) extends FileProvider {
  require(trackName != null, "RodBind trackName cannot be null")
  require(trackType != null, "RodBind trackType cannot be null")
  require(file != null, "RodBind file cannot be null")
}

/**
 * Used to provide -B rodBind arguments to the GATK.
 */
object RodBind {
  /**
   * Formats the rod binding on the command line.
   * Used for optional and repeat.
   * @param cmdLineParam command line parameter, ex: -B
   * @param prefix unused
   * @param value RodBind to add.
   * @param suffix unused
   * @return The command line addition.
   */
  def formatCommandLine(cmdLineParam: String)(prefix: String, value: Any, suffix: String) = {
    value match {
      case rodBind: RodBind =>
        " %s:%s,%s %s".format(cmdLineParam, rodBind.trackName, rodBind.trackType, rodBind.file)
    }
  }
}
