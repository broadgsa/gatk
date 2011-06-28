package org.broadinstitute.sting.queue.util

import org.apache.commons.io.FileUtils
import java.io.{FileReader, File}
import org.broadinstitute.sting.utils.exceptions.UserException

/**
 * A collection of utilities for modifying java.io.
 */
object IOUtils extends Logging {
  /**
   * Checks if the temp directory has been setup and throws an exception if they user hasn't set it correctly.
   * @param tempDir Temporary directory.
   */
  def checkTempDir(tempDir: File) = {
    val tempDirPath = tempDir.getAbsolutePath
    // Keeps the user from leaving the temp directory as the default, and on Macs from having pluses
    // in the path which can cause problems with the Google Reflections library.
    // see also: http://benjchristensen.com/2009/09/22/mac-osx-10-6-java-java-io-tmpdir/
    if (tempDirPath.startsWith("/var/folders/") || (tempDirPath == "/tmp") || (tempDirPath == "/tmp/"))
      throw new UserException.BadTmpDir("java.io.tmpdir must be explicitly set")
    if (!tempDir.exists && !tempDir.mkdirs)
      throw new UserException.BadTmpDir("Could not create directory: " + tempDir.getAbsolutePath())
  }

  /**
   * Creates a temp directory with the prefix and optional suffix.
   * @param prefix Prefix for the directory name.
   * @param suffix Optional suffix for the directory name.
   * @param tempDirParent Parent directory for the temp directory.
   * @return The created temporary directory.
   */
  def tempDir(prefix: String, suffix: String = "", tempDirParent: File) = {
    if (!tempDirParent.exists && !tempDirParent.mkdirs)
       throw new UserException.BadTmpDir("Could not create temp directory: " + tempDirParent)
    val temp = File.createTempFile(prefix + "-", suffix, tempDirParent)
    if (!temp.delete)
      throw new UserException.BadTmpDir("Could not delete sub file: " + temp.getAbsolutePath())
    if (!temp.mkdir)
      throw new UserException.BadTmpDir("Could not create sub directory: " + temp.getAbsolutePath())
    absolute(temp)
  }

  /**
   * Writes content into a file.
   * @param file File to write to.
   * @param content Content to write.
   */
  def writeContents(file: File, content: String) =  FileUtils.writeStringToFile(file, content)

  /**
   * Reads content of a file into a string.
   * Only for use on really small files!
   * @param file File to read to.
   * @return content Content of the file.
   */
  def readContents(file: File) = FileUtils.readFileToString(file)

  /**
   * Writes content to a temp file and returns the path to the temporary file.
   * @param content to write.
   * @param prefix Prefix for the temp file.
   * @param suffix Suffix for the temp file.
   * @param directory Directory for the temp file.
   * @return the path to the temp file.
   */
  def writeTempFile(content: String, prefix: String, suffix: String, directory: File) = {
    val tempFile = absolute(File.createTempFile(prefix, suffix, directory))
    writeContents(tempFile, content)
    tempFile
  }

  /**
   * Waits for NFS to propagate a file creation, imposing a timeout.
   *
   * Based on Apache Commons IO FileUtils.waitFor()
   *
   * @param file The file to wait for.
   * @param seconds The maximum time in seconds to wait.
   * @return true if the file exists
   */
  def waitFor(file: File, seconds: Int): Boolean = waitFor(List(file), seconds).isEmpty

  /**
   * Waits for NFS to propagate a file creation, imposing a timeout.
   *
   * Based on Apache Commons IO FileUtils.waitFor()
   *
   * @param files The list of files to wait for.
   * @param seconds The maximum time in seconds to wait.
   * @return Files that still do not exists at the end of the timeout, or a empty list if all files exists.
   */
  def waitFor[T <: Traversable[File]](files: T, seconds: Int): Traversable[File] = {
      var timeout = 0;
      var tick = 0;
      var missingFiles = files.filterNot(_.exists)
      while (!missingFiles.isEmpty && timeout <= seconds) {
        if (tick >= 10) {
          tick = 0;
          timeout += 1
        }
        tick += 1
        try {
            Thread.sleep(100)
        } catch {
          case ignore: InterruptedException =>
        }
        missingFiles = missingFiles.filterNot(_.exists)
      }
      missingFiles
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
   * Returns the sub path rooted at the parent.
   * @param parent The parent directory.
   * @param path The sub path to append to the parent, if the path is not absolute.
   * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
   */
  def absolute(parent: File, path: String): File =
    absolute(parent, new File(path))

  /**
   * Returns the sub path rooted at the parent.
   * @param parent The parent directory.
   * @param file The sub path to append to the parent, if the path is not absolute.
   * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
   */
  def absolute(parent: File, file: File): File = {
    if (file.isAbsolute)
      absolute(file)
    else
      absolute(new File(parent, file.getPath))
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

  /**
   * Tries to delete a file.  Emits a warning if the file was unable to be deleted.
   * @param file File to delete.
   * @return true if the file was deleted.
   */
  def tryDelete(file: File) = {
    val deleted = FileUtils.deleteQuietly(file)
    if (deleted)
      logger.debug("Deleted " + file)
    else if (file.exists)
      logger.warn("Unable to delete " + file)
    deleted
  }
}
