package org.broadinstitute.sting.queue.function

import java.io.File

trait IntervalFunction extends InputOutputFunction {
  var referenceFile: File
  var intervals: File
}
