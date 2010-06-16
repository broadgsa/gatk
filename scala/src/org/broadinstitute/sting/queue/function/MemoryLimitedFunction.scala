package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.util.{Input, Optional}

trait MemoryLimitedFunction {
  @Input
  @Optional
  var memoryLimit: Option[Int] = None
}
