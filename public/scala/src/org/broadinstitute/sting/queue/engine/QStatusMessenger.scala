package org.broadinstitute.sting.queue.engine

/**
 * Plugin to sends QStatus messages
 */
trait QStatusMessenger {
  def started()
  def done()
  def exit(message: String)
}
