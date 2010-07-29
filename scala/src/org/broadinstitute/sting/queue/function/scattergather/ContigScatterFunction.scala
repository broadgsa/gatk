package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.Input
import org.broadinstitute.sting.queue.function.IntervalFunction

class ContigScatterFunction extends ScatterFunction {
  type ScatterType = File

  @Input(doc="Reference file to scatter")
  var referenceFile: File = _

  override def setOriginalFunction(originalFunction: ScatterGatherableFunction) = {
    val command = originalFunction.asInstanceOf[IntervalFunction]
    referenceFile = command.referenceFile
    super.setOriginalFunction(originalFunction)
  }

  // TODO: Use the reference file for "all"
  def commandLine = "splitIntervalsByContig.py %s%s".format(originalInput, repeat(" ", scatterParts))
}
