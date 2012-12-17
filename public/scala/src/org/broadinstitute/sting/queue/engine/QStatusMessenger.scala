package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.util.RemoteFile

/**
 * Plugin to sends QStatus messages
 */
trait QStatusMessenger {
  def started()
  def done(inputs: Seq[_], outputs: Seq[_])
  def exit(message: String)

  def started(job: String)
  def done(job: String)
  def exit(job: String, message: String)
}
