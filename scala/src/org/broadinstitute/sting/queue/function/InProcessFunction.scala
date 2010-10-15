package org.broadinstitute.sting.queue.function

/**
 * Runs a function in process.
 */
trait InProcessFunction extends QFunction {
  def run()
  def description = (List(this.getClass.getSimpleName) ++ this.outputs.filter(file => useStatusOutput(file)).map(_.getAbsolutePath)).mkString(" ")
}
