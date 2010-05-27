package org.broadinstitute.sting.queue.engine

import org.apache.commons.lang.builder.{EqualsBuilder, HashCodeBuilder}
import java.io.File
import org.apache.commons.lang.StringUtils
import org.broadinstitute.sting.queue.QException

/**
 * Represents a file extension along with several tags.
 * TODO: Use the tags to map rules between wildcards, ex: *.vcf -> *.eval
 */
class QFile(val fileType: String, val parts: String*) {
  val extension = (List(parts:_*) ::: List(fileType)).mkString(".")
  override def toString = extension
  override def equals(p1: Any) = EqualsBuilder.reflectionEquals(this, p1)
  override def hashCode = HashCodeBuilder.reflectionHashCode(this)

  def matchesFile(path: String): Boolean = matchesFile(new File(path))
  def matchesFile(file: File): Boolean = file.getCanonicalPath.endsWith(extension)
  def baseName(path: String): String = baseName(new File(path))
  def baseName(file: File): String = StringUtils.removeEnd(file.getAbsolutePath, extension)
  def fullName(baseName: String) = baseName + extension
}

object QFile {
  def getFiles(files: Any) : List[QFile] = {
    files match {
      case null => List.empty[QFile]
      case Nil => List.empty[QFile]
      case path: String => List(new QFile(path))
      case file: QFile => List(file)
      // Any List or Tuple add the members to this list
      case product: Product => {
        var list = List.empty[QFile]
        for (fileList <- product.productIterator.toList.map(getFiles(_))) {
          list :::= fileList
        }
        list
      }
      case x => throw new QException("Unknown file type: " + x)
    }
  }
}
