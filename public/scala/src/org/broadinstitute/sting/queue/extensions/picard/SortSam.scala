package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File
import org.broadinstitute.sting.queue.QScript._

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class SortSam extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "SortSam"
  javaMainClass = "net.sf.picard.sam.SortSam"

  @Input(doc="The input SAM or BAM files to sort.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: List[File] = Nil

  @Output(doc="The sorted BAM or SAM output file.", shortName = "output", fullName = "output_bam_file", required = true)
  var output: File = _

  @Output(doc="The output bam index", shortName = "out_index", fullName = "output_bam_index_file", required = false)
  var outputIndex: File = _

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }



  override def inputBams = input
  override def outputBam = output
  this.createIndex = Some(true)
  override def commandLine = super.commandLine
}