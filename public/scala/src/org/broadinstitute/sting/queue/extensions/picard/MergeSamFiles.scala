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
class MergeSamFiles extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "MergeSamFiles"
  javaMainClass = "net.sf.picard.sam.MergeSamFiles"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: List[File] = Nil

  @Output(doc="The output merged BAM file", shortName = "output", fullName = "output_bam_file", required = true)
  var output: File = _

  @Output(doc="The output bam index", shortName = "out_index", fullName = "output_bam_index_file", required = false)
  var outputIndex: File = _

  @Argument(doc="Merge the seqeunce dictionaries Default value: false. This option can be set to 'null' to clear the default value.", shortName = "merge_dict", fullName = "merge_sequence_dictionaries", required = false)
  var MERGE_SEQUENCE_DICTIONARIES: Boolean = false

  @Argument(doc = "Option to enable a simple two-thread producer consumer version of the merge algorithm that uses one thread to read and merge the records from the input files and another thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file. ", shortName = "thread", fullName ="use_threading", required=false)
  var USE_THREADING: Boolean = false;

  @Argument(doc = "Comments to include in the merged output file's header.", shortName = "com", fullName = "comments", required = false)
  var COMMENT: String = ""

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }


  override def inputBams = input
  override def outputBam = output
  this.createIndex = Some(true)
  override def commandLine = super.commandLine +
       conditionalParameter(MERGE_SEQUENCE_DICTIONARIES, " MERGE_SEQUENCE_DICTIONARIES=true") +
       conditionalParameter(USE_THREADING, " USE_THREADING=true") +
       conditionalParameter(COMMENT != null && !COMMENT.isEmpty, " COMMENT=" + COMMENT)


}