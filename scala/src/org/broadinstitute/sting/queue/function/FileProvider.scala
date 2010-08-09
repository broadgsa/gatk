package org.broadinstitute.sting.queue.function

import java.io.File

/**
 * An trait for @Input or @Output CommandLineFunction fields that are not files, but have a File that can be get/set.
 */
trait FileProvider {
  /** Gets/Sets the file. */
  var file: File
}
