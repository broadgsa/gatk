package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class MarkDuplicates extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "MarkDuplicates"
  javaMainClass = "net.sf.picard.sam.MarkDuplicates"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: List[File] = Nil

  @Output(doc="The output file to write marked records to", shortName = "output", fullName = "output_bam_file", required = true)
  var output: File = _

  @Output(doc="The output bam index", shortName = "out_index", fullName = "output_bam_index_file", required = false)
  var outputIndex: File = _

  @Output(doc="File to write duplication metrics to", shortName = "out_metrics", fullName = "output_metrics_file", required = false)
  var metrics: File = new File(output + ".metrics")

  @Argument(doc="If true do not write duplicates to the output file instead of writing them with appropriate flags set.", shortName = "remdup", fullName = "remove_duplicates", required = false)
  var REMOVE_DUPLICATES: Boolean = false

  @Argument(doc = "Maximum number of file handles to keep open when spilling read ends to disk.  Set this number a little lower than the per-process maximum number of file that may be open.  This number can be found by executing the 'ulimit -n' command on a Unix system.", shortName = "max_file_handles", fullName ="max_file_handles_for_read_ends_maps", required=false)
  var MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: Int = -1;

  @Argument(doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections.  If you are running out of memory, try reducing this number.", shortName = "sorting_ratio", fullName = "sorting_collection_size_ratio", required = false)
  var SORTING_COLLECTION_SIZE_RATIO: Double = -1

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }


  override def inputBams = input
  override def outputBam = output
  this.sortOrder = null
  this.createIndex = Some(true)
  override def commandLine = super.commandLine +
       " M=" + metrics +
       conditionalParameter(REMOVE_DUPLICATES, " REMOVE_DUPLICATES=true") +
       conditionalParameter(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP > 0, " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=" + MAX_FILE_HANDLES_FOR_READ_ENDS_MAP.toString) +
       conditionalParameter(SORTING_COLLECTION_SIZE_RATIO > 0, " SORTING_COLLECTION_SIZE_RATIO=" + SORTING_COLLECTION_SIZE_RATIO.toString)


}