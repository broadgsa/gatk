package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.commandline.{Argument, Output, Input}

/**
 * Indexes a BAM file.
 * By default uses samtools index.
 * The syntax of the script must be:
 * <bamIndexScript> <bam_file> <bam_index_file>
 */
class BamIndexFunction extends CommandLineFunction {
  @Argument(doc="BAM file script")
  var bamIndexScript: String = "samtools index"

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

  def commandLine = "%s %s %s".format(bamIndexScript, bamFile, bamFileIndex)

  override def dotString = "Index: %s".format(bamFile.getName)
}
