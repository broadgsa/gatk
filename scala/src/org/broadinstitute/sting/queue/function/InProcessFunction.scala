package org.broadinstitute.sting.queue.function

import java.io.File

/**
 * Runs a function in process.
 */
trait InProcessFunction extends QFunction {
  def run()
  def useStatusOutput(file: File) = true
  def description = (List(this.getClass.getSimpleName) ++ this.outputs.map(_.getAbsolutePath)).mkString(" ")
}
