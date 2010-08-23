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
    subDir(dir, new File(path))

  /**
   * Returns the sub path rooted at the parent.
   * If the sub path is already absolute, returns the sub path.
   * If the parent is the current directory, returns the sub path.
   * If the sub path is the current directory, returns the parent.
   * Else returns new File(parent, subPath)
   * @param parent The parent directory
   * @param file The sub path to append to the parent, if the path is not absolute.
   * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
   */
  def subDir(parent: File, file: File): File = {
    val parentAbs = absolute(parent)
    val fileAbs = absolute(file)
    val currentAbs = absolute(CURRENT_DIR)
    if (parentAbs == currentAbs && fileAbs == currentAbs)
      absolute(CURRENT_DIR.getCanonicalFile)
    else if (parentAbs == currentAbs || file.isAbsolute)
      fileAbs
    else if (fileAbs == currentAbs)
      parentAbs
    else
      absolute(new File(parentAbs, file.getPath))
  }

  /**
   * Resets the parent of the file to the directory.
   * @param dir New parent directory.
   * @param file Path to the file to be re-rooted.
   * @return Absolute path to the new file.
   */
  def resetParent(dir: File, file: File) = absolute(subDir(dir, file.getName))

  /**
   * Creates a scatterGatherTempDir directory with the prefix and optional suffix.
   * @param prefix Prefix for the directory name.
   * @param suffix Optional suffix for the directory name.  Defaults to "".
   * @return The created temporary directory.
   * @throws IOException if the directory could not be created.
   */
  def tempDir(prefix: String, suffix: String = "") = {
    val tempDirParent = new File(System.getProperty("java.io.tmpdir"))
    if (!tempDirParent.exists && !tempDirParent.mkdirs)
       throw new IOException("Could not create temp directory: " + tempDirParent)
    val temp = File.createTempFile(prefix + "-", suffix)
    if (!temp.delete)
      throw new IOException("Could not delete sub file: " + temp.getAbsolutePath())
    if (!temp.mkdir)
      throw new IOException("Could not create sub directory: " + temp.getAbsolutePath())
    absolute(temp)
  }

  /**
   * Returns the directory at the number of levels deep.
   * For example 2 levels of /path/to/dir will return /path/to
   * @param dir Directory path.
   * @param level how many levels deep from the root.
   * @return The path to the parent directory that is level-levels deep.
   */
  def dirLevel(dir: File, level: Int): File = {
    var directories = List.empty[File]
    var parentDir = absolute(dir)
    while (parentDir != null) {
      directories +:= parentDir
      parentDir = parentDir.getParentFile
    }
    if (directories.size <= level)
      directories.last
    else
      directories(level)
  }

  def absolute(file: File) = {
    var fileAbs = file.getAbsoluteFile
    var names = List.empty[String]
    while (fileAbs != null) {
      val name = fileAbs.getName
      fileAbs = fileAbs.getParentFile
      
      if (name == ".") {
        /* skip */

        /* TODO: What do we do for ".."?
      } else if (name == "..") {

        CentOS tcsh says use getCanonicalFile:
        ~ $ mkdir -p test1/test2
        ~ $ ln -s test1/test2 test3
        ~ $ cd test3/..
        ~/test1 $

        Mac bash says keep going with getAbsoluteFile:
        ~ $ mkdir -p test1/test2
        ~ $ ln -s test1/test2 test3
        ~ $ cd test3/..
        ~ $

        For now, leave it and let the shell figure it out.
        */
      } else {
        names +:= name
      }
    }

    new File(names.mkString("/", "/", ""))
  }
}
