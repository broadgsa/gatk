package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import java.io.File

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 10/9/12
 * Time: 5:59 PM
 * To change this template use File | Settings | File Templates.
 */
class CalculateHsMetrics extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "CalculateHsMetrics"
  javaMainClass = "net.sf.picard.sam.CalculateHsMetrics"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

  @Output(doc="The output file to write statistics to", shortName = "output", fullName = "output_file", required = true)
  var output: File = _

  @Argument(doc="Interval list with targets", shortName = "targets", fullName = "target_list", required = true)
  var targets: File = _

  @Argument(doc="Interval list with baits", shortName = "baits", fullName = "bait_list", required = true)
  var baits: File = _

  @Argument(doc="Reference file", shortName = "reference", fullName = "reference", required = true)
  var reference: File = _
 /*
  @Argument(doc = "Maximum number of file handles to keep open when spilling read ends to disk.  Set this number a little lower than the per-process maximum number of file that may be open.  This number can be found by executing the 'ulimit -n' command on a Unix system.", shortName = "max_file_handles", fullName ="max_file_handles_for_read_ends_maps", required=false)
  var MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: Int = -1;

  @Argument(doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections.  If you are running out of memory, try reducing this number.", shortName = "sorting_ratio", fullName = "sorting_collection_size_ratio", required = false)
  var SORTING_COLLECTION_SIZE_RATIO: Double = -1
   */
  override def freezeFieldValues() {
    super.freezeFieldValues()
//    if (outputIndex == null && output != null)
  //    outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }

  val level = "SAMPLE"

  override def inputBams = input
  override def outputBam = output
  //this.sortOrder = null
  //this.createIndex = Some(true)
  override def commandLine = super.commandLine +
    required("BAIT_INTERVALS=" + baits) +
    required("TARGET_INTERVALS=" + targets) +
    required("REFERENCE_SEQUENCE=" + reference) +
    optional("METRIC_ACCUMULATION_LEVEL="+level)/*+
    conditional(REMOVE_DUPLICATES, "REMOVE_DUPLICATES=true") +
    conditional(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP > 0, "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=" + MAX_FILE_HANDLES_FOR_READ_ENDS_MAP.toString) +
    conditional(SORTING_COLLECTION_SIZE_RATIO > 0, "SORTING_COLLECTION_SIZE_RATIO=" + SORTING_COLLECTION_SIZE_RATIO.toString)    */


}
