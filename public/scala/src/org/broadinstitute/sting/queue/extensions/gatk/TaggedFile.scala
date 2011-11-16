package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.utils.io.FileExtension

/**
 * Used to provide tagged -I input_file arguments to the GATK.
 */
class TaggedFile(path: String, val tag: String) extends File(path) with FileExtension {
  def this(file: File, tag: String) =
    this(file.getPath, tag)
  def withPath(path: String) = new TaggedFile(path, tag)
}

/**
 * Used to provide -I input_file arguments to the GATK.
 */
object TaggedFile {
  def apply(path: String, tag: String) = new TaggedFile(path, tag)
  def apply(file: File, tag: String) = new TaggedFile(file, tag)

  /**
   * Formats the rod binding on the command line.
   * Used for optional and repeat.
   * @param cmdLineParam command line parameter, ex: -I
   * @param prefix unused
   * @param value TaggedFile to add.
   * @param suffix unused
   * @return The command line addition.
   */
  def formatCommandLine(cmdLineParam: String)(prefix: String, value: Any, suffix: String) = {
    value match {
      case taggedFile: TaggedFile if (taggedFile.tag != null) =>
          " %s:%s %s".format(cmdLineParam, taggedFile.tag, taggedFile.getPath)
      case file: File =>
          " %s %s".format(cmdLineParam, file.getPath)
    }
  }
}
