package org.broadinstitute.sting.queue.extensions.samtools

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.commandline.{Argument, Output, Input}

/**
 * Indexes a BAM file using samtools.
 */
class SamtoolsIndexFunction extends CommandLineFunction {
  @Argument(doc="samtools path")
  var samtools: String = "samtools"

  @Input(doc="BAM file to index")
  var bamFile: File = _

  @Output(doc="BAM file index to output", required=false)
  var bamFileIndex: File = _

  /**
   * Sets the bam file index to the bam file name + ".bai".
   */
  override def freezeFieldValues = {
    super.freezeFieldValues
    if (bamFileIndex == null && bamFile != null)
      bamFileIndex = new File(bamFile.getPath + ".bai")
  }

  def commandLine = "%s index %s %s".format(samtools, bamFile, bamFileIndex)

  override def dotString = "Index: %s".format(bamFile.getName)
}
