package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.FileProvider

/**
 * Used to provide -B rodBinding arguments to the GATK.
 */
case class RodBind(var trackName: String, var trackType: String, var file: File) extends FileProvider {
  require(trackName != null, "RodBind trackName cannot be null")
  require(trackType != null, "RodBind trackType cannot be null")
  require(file != null, "RodBind file cannot be null")
  override def toString = "%s,%s,%s".format(trackName, trackType, file)
}
