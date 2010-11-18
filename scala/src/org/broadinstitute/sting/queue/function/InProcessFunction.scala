package org.broadinstitute.sting.queue.function

/**
 * Runs a function in process.
 */
trait InProcessFunction extends QFunction {
  def run()
  def description = (List(this.getClass.getSimpleName) ++ this.outputs.filterNot(file => isLogFile(file)).map(_.getAbsolutePath)).mkString(" ")
}
