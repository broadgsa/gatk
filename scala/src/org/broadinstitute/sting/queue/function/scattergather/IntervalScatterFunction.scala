package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.util.Input
import org.broadinstitute.sting.queue.function.IntervalFunction

class IntervalScatterFunction extends ScatterFunction {
  type ScatterType = File

  @Input
  var referenceFile: File = _

  override def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {
    val command = originalFunction.asInstanceOf[IntervalFunction]
    referenceFile = command.referenceFile
    super.setOriginalFunction(originalFunction)
  }

  // TODO: Use the reference file for "all"
  def commandLine = "splitIntervals.sh %s%s".format(originalInput, repeat(" ", scatterParts))
}
