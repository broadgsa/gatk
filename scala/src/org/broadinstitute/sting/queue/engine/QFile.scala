package org.broadinstitute.sting.queue.engine

import org.apache.commons.lang.builder.{EqualsBuilder, HashCodeBuilder}
import java.io.File
import org.apache.commons.lang.StringUtils

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
  def baseName(file: File): String = StringUtils.removeEnd(file.getCanonicalPath, extension)
  def fullName(baseName: String) = baseName + extension
}
