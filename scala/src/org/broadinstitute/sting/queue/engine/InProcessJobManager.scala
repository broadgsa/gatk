package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.InProcessFunction

class InProcessJobManager extends JobManager[InProcessFunction, InProcessRunner] {
  def runnerType = classOf[InProcessRunner]
  def functionType = classOf[InProcessFunction]
  def create(function: InProcessFunction) = new InProcessRunner(function)
}
