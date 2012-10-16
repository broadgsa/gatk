package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.util.RemoteFile

/**
 * Plugin to sends QStatus messages
 */
trait QStatusMessenger {
  def started()
  def done(files: Seq[Map[ArgumentSource, Seq[RemoteFile]]])
  def exit(message: String)
}
