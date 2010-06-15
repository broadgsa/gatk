package org.broadinstitute.sting.queue.function

trait IntervalFunction {
  type Intervals = String
  var intervals: Intervals
}