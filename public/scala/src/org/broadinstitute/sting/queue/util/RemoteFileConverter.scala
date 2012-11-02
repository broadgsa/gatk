package org.broadinstitute.sting.queue.util

import java.io.File

trait RemoteFileConverter {
  type RemoteFileType <: RemoteFile

  /**
   * If this remote file creator is capable of converting to a remote file.
   * @return true if ready to convert
   */
  def convertToRemoteEnabled: Boolean

  /**
   * Converts to a remote file
   * @param file The original file
   * @param name A "name" to use for the remote file
   * @return The new version of this remote file.
   */
  def convertToRemote(file: File, name: String): RemoteFileType
}
