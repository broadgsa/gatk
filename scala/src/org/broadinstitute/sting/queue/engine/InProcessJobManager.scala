package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.InProcessFunction

class InProcessJobManager extends JobManager[InProcessFunction, InProcessRunner] {
  def create(function: InProcessFunction) = new InProcessRunner(function)
}
