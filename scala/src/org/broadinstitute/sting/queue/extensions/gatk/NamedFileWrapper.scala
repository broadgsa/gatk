package org.broadinstitute.sting.queue.extensions.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.FileProvider

class NamedFileWrapper(private val file: File) {
  def toNamedFile = new NamedFile(file)
  def toNamedFile(name: String) = new NamedFile(file, name)
}
