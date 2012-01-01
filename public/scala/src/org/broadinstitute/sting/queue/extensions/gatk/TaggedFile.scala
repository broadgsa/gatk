package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.utils.io.FileExtension
import org.broadinstitute.sting.queue.util.ShellUtils

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

  def formatCommandLineParameter(cmdLineParam: String, value: Any) = {
    value match {
      case taggedFile: TaggedFile if (taggedFile.tag != null) =>
        "%s:%s".format(cmdLineParam, taggedFile.tag)
      case file: File =>
        cmdLineParam
      case x =>
        ""
    }
  }
}
