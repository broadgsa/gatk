package org.broadinstitute.sting.queue.util

import java.io.{IOException, File}

object IOUtils {
  val CURRENT_DIR = new File(".")

  def sub(parent: File, subPath: String) = {
    val file = new File(subPath)
    if (parent == CURRENT_DIR && file == CURRENT_DIR)
      CURRENT_DIR.getCanonicalFile
    else if (parent == CURRENT_DIR || file.isAbsolute)
      file
    else if (file == CURRENT_DIR)
      parent
    else
      new File(parent, subPath)
  }

  def temp(prefix: String, suffix: String = "") = {
    val tempDir = File.createTempFile(prefix + "-", suffix)
    if(!tempDir.delete)
      throw new IOException("Could not delete sub file: " + tempDir.getAbsolutePath())
    if(!tempDir.mkdir)
      throw new IOException("Could not create sub directory: " + tempDir.getAbsolutePath())
    tempDir
  }

  def reset(dir: File, file: File) = sub(dir, file.getName).getAbsoluteFile
  def absolute(dir: File, file: File) = sub(dir, file.getPath).getAbsoluteFile
}
