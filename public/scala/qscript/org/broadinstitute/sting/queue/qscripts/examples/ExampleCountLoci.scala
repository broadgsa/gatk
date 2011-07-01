package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An introductory pipeline with integration tests testing the MD5 of the @Output.
 */
class ExampleCountLoci extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = null

  @Input(doc="One or more bam files.", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="Intervals to traverse.", shortName="L", required=false)
  var intervals: List[String] = Nil

  @Output
  var out: File = _

  def script = {
    val countLoci = new CountLoci
    countLoci.reference_sequence = referenceFile
    countLoci.input_file = bamFiles
    countLoci.intervalsString = intervals
    countLoci.out = out
    countLoci.memoryLimit = 1
    add(countLoci)
  }
}
