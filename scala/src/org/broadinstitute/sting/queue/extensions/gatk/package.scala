package org.broadinstitute.sting.queue.extensions

import java.io.File
import org.broadinstitute.sting.queue.extensions.gatk.NamedFile
import org.broadinstitute.sting.queue.extensions.gatk.NamedFileWrapper

package object gatk {
  implicit def fileToNamedFileWrapper(file: File) = new NamedFileWrapper(file)
  // TODO: Get the syntax right so that the implicits kick in for a generic type, ex: Travesable[File], Traversable[_ <: File], etc.
  // but need to return the same outter type, so T <: Traversable[File] : T[NamedFile], T <: Traversable[_ <: File]: T[NamedFile], etc.
  implicit def filesToNamedFilesWrapper(files: List[File]) = files.map(file => if (file == null) null else new NamedFile(file))
  implicit def filesToNamedFilesWrapper(files: Set[File]) = files.map(file => if (file == null) null else new NamedFile(file))
}
