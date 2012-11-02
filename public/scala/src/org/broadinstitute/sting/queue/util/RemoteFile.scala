package org.broadinstitute.sting.queue.util

import java.io.File
import org.broadinstitute.sting.utils.io.FileExtension

/**
 * An extension of java.io.File that can be pulled from or pushed to a remote location.
 */
trait RemoteFile extends File with FileExtension {
  def pullToLocal()
  def pushToRemote()
  def deleteRemote()
  def remoteDescription: String
}
