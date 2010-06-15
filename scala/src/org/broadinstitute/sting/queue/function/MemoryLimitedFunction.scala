package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.util.{Input, Optional, ClassType}

trait MemoryLimitedFunction {
  @Input
  @Optional
  @ClassType(classOf[Int])
  var memoryLimit: Option[Int] = None
}
