package org.broadinstitute.sting.queue.function

import java.io.File

/**
 * An trait for @Input or @Output CommandLineFunction fields that are extensions of files.
 */
trait FileExtension extends File {
  /**
   * Returns a clone of the FileExtension with the new path.
   * @param newPath new path for the clone of this FileExtension
   * @return a clone of the FileExtension with the new path.
   */
  def withPath(newPath: String): File
}
