package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.FileExtension
import java.lang.String

/**
 * Used to provide -B rodBind arguments to the GATK.
 */
class RodBind(var trackName: String, var trackType: String, path: String) extends File(path) with FileExtension {
  def this(trackName: String, trackType: String, file: File) =
    this(trackName, trackType, file.getPath)
  require(trackName != null, "RodBind trackName cannot be null")
  require(trackType != null, "RodBind trackType cannot be null")
  def withPath(newPath: String) = new RodBind(trackName, trackType, newPath)
}

/**
 * Used to provide -B rodBind arguments to the GATK.
 */
object RodBind {
  def apply(trackName: String, trackType: String, path: String) = new RodBind(trackName, trackType, path)
  def apply(trackName: String, trackType: String, file: File) = new RodBind(trackName, trackType, file)

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
        " %s:%s,%s %s".format(cmdLineParam, rodBind.trackName, rodBind.trackType, rodBind.getPath)
    }
  }
}
