package org.broadinstitute.sting.queue.util

import java.io.{IOException, File}

/**
 * A collection of utilities for modifying java.io.
 */
object IOUtils {
  /** The current directory "." */
  val CURRENT_DIR = new File(".")


  /**
   * Returns the sub path rooted at the parent.
   * If the sub path is already absolute, returns the sub path.
   * If the parent is the current directory, returns the sub path.
   * If the sub bath is the current directory, returns the parent.
   * Else returns new File(parent, subPath)
   * @param parent The parent directory
   * @param path The sub path to append to the parent, if the path is not absolute.
   * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
   */
  def subDir(dir: File, path: String): File =
    subDir(dir.getAbsoluteFile, new File(path))

  /**
   * Returns the sub path rooted at the parent.
   * If the sub path is already absolute, returns the sub path.
   * If the parent is the current directory, returns the sub path.
   * If the sub bath is the current directory, returns the parent.
   * Else returns new File(parent, subPath)
   * @param parent The parent directory
   * @param file The sub path to append to the parent, if the path is not absolute.
   * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
   */
  def subDir(parent: File, file: File): File = {
    if (parent == CURRENT_DIR && file == CURRENT_DIR)
      CURRENT_DIR.getCanonicalFile.getAbsoluteFile
    else if (parent == CURRENT_DIR || file.isAbsolute)
      file.getAbsoluteFile
    else if (file == CURRENT_DIR)
      parent.getAbsoluteFile
    else
      new File(parent, file.getPath).getAbsoluteFile
  }

  /**
   * Resets the parent of the file to the directory.
   * @param dir New parent directory.
   * @param file Path to the file to be re-rooted.
   * @return Absolute path to the new file.
   */
  def resetParent(dir: File, file: File) = subDir(dir.getAbsoluteFile, file.getName).getAbsoluteFile

  /**
   * Creates a scatterGatherTempDir directory with the prefix and optional suffix.
   * @param prefix Prefix for the directory name.
   * @param suffix Optional suffix for the directory name.  Defaults to "".
   * @return The created temporary directory.
   * @throws IOException if the directory could not be created.
   */
  def tempDir(prefix: String, suffix: String = "") = {
    val temp = File.createTempFile(prefix + "-", suffix)
    if(!temp.delete)
      throw new IOException("Could not delete sub file: " + temp.getAbsolutePath())
    if(!temp.mkdir)
      throw new IOException("Could not create sub directory: " + temp.getAbsolutePath())
    temp
  }
}
