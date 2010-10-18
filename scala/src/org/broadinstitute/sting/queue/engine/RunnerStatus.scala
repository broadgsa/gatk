package org.broadinstitute.sting.queue.engine

object RunnerStatus extends Enumeration {
  val PENDING = Value("pending")
  val RUNNING = Value("running")
  val FAILED = Value("failed")
  val DONE = Value("done")
  val SKIPPED = Value("skipped")
}
