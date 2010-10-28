package org.broadinstitute.sting.queue.util

import org.apache.commons.io.FileUtils
import java.io.{FileReader, File}
import org.broadinstitute.sting.utils.exceptions.UserException

/**
 * A collection of utilities for modifying java.io.
 */
object IOUtils {
  /** The current directory "." */
  val CURRENT_DIR = new File(".")

  val CURRENT_DIR_ABS = absolute(CURRENT_DIR)

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
  def subDir(parent: File, path: String): File =
    subDir(parent, new File(path))

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
    if (parentAbs == CURRENT_DIR_ABS && fileAbs == CURRENT_DIR_ABS)
      CURRENT_DIR_ABS
    else if (parentAbs == CURRENT_DIR_ABS || file.isAbsolute)
      fileAbs
    else if (fileAbs == CURRENT_DIR_ABS)
      parentAbs
    else
      absolute(new File(parentAbs, file.getPath))
  }

  def checkTempDir = {
    val javaTemp = System.getProperty("java.io.tmpdir")
    // Keeps the user from leaving the temp directory as the default, and on Macs from having pluses
    // in the path which can cause problems with the Google Reflections library.
    // see also: http://benjchristensen.com/2009/09/22/mac-osx-10-6-java-java-io-tmpdir/
    if (javaTemp.startsWith("/var/folders/") || (javaTemp == "/tmp") || (javaTemp == "/tmp/"))
      throw new UserException.BadTmpDir("java.io.tmpdir must be explicitly set")
  }

  /**
   * Returns the temp directory as defined by java.
   * @return the temp directory as defined by java.
   */
  def javaTempDir() = {
    val tempDir = new File(System.getProperty("java.io.tmpdir"))
    if (!tempDir.exists && !tempDir.mkdirs)
      throw new UserException.BadTmpDir("Could not create directory: " + tempDir.getAbsolutePath())
    absolute(tempDir)
  }

  /**
   * Creates a temp directory with the prefix and optional suffix.
   * @param prefix Prefix for the directory name.
   * @param suffix Optional suffix for the directory name.  Defaults to "".
   * @return The created temporary directory.
   */
  def tempDir(prefix: String, suffix: String = "") = {
    val tempDirParent = javaTempDir()
    if (!tempDirParent.exists && !tempDirParent.mkdirs)
       throw new UserException.BadTmpDir("Could not create temp directory: " + tempDirParent)
    val temp = File.createTempFile(prefix + "-", suffix)
    if (!temp.delete)
      throw new UserException.BadTmpDir("Could not delete sub file: " + temp.getAbsolutePath())
    if (!temp.mkdir)
      throw new UserException.BadTmpDir("Could not create sub directory: " + temp.getAbsolutePath())
    absolute(temp)
  }

  def writeContents(file: File, content: String) =  FileUtils.writeStringToFile(file, content)

  def   writeTempFile(content: String, prefix: String, suffix: String = "", directory: File = null) = {
    val tempFile = absolute(File.createTempFile(prefix, suffix, directory))
    writeContents(tempFile, content)
    tempFile
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

  /**
   * A mix of getCanonicalFile and getAbsoluteFile that returns the
   * absolute path to the file without deferencing symbolic links.
   * @param file the file.
   * @return the absolute path to the file.
   */
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

  /**
   * Returns the last lines of the file.
   * NOTE: This is only safe to run on smaller files!
   * @param file File to read.
   * @param count Maximum number of lines to return.
   * @return The last count lines from file.
   */
  def tail(file: File, count: Int) = {
    var tailLines = List.empty[String]
    var reader = new FileReader(file)
    try {
      val iterator = org.apache.commons.io.IOUtils.lineIterator(reader)
      var lineCount = 0
      while (iterator.hasNext) {
        val line = iterator.nextLine
        lineCount += 1
        if (lineCount > count)
          tailLines = tailLines.tail
        tailLines :+= line
      }
    } finally {
      org.apache.commons.io.IOUtils.closeQuietly(reader)
    }
    tailLines
  }
}
