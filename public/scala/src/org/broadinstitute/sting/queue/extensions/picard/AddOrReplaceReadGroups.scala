package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class AddOrReplaceReadGroups extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "AddOrReplaceReadGroups"
  javaMainClass = "net.sf.picard.sam.AddOrReplaceReadGroups"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: List[File] = Nil

  @Output(doc="The output BAM file with the modified/added read groups", shortName = "output", fullName = "output_bam_file", required = true)
  var output: File = _

  @Output(doc="The output bam index", shortName = "out_index", fullName = "output_bam_index_file", required = false)
  var outputIndex: File = _

  @Argument(doc="Read group ID", shortName = "id", fullName = "read_group_id", required = true)
  var RGID: String = _

  @Argument(doc = "Read group library", shortName = "lb", fullName = "read_group_library", required = true)
  var RGLB: String = _

  @Argument(doc = "Read group platform (e.g. illumina, solid)", shortName = "pl", fullName ="read_group_platform", required=true)
  var RGPL: String = _

  @Argument(doc = "Read group platform unit (e.g. run barcode)", shortName = "pu", fullName = "read_group_platform_unit", required = true)
  var RGPU: String = _

  @Argument(doc = "Read group sample name", shortName = "sm", fullName = "read_group_sample_name", required = true)
  var RGSM: String = _

  @Argument(doc = "Read group center name", shortName = "cn", fullName = "read_group_center_name", required = false)
  var RGCN: String = ""

  @Argument(doc = "Read group description", shortName = "ds", fullName = "read_group_description", required = false)
  var RGDS: String = ""

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }


  override def inputBams = input
  override def outputBam = output
  this.createIndex = Some(true)
  override def commandLine = super.commandLine +
       " RGID=" + RGID +
       " RGLB=" + RGLB +
       " RGPL=" + RGPL +
       " RGPU=" + RGPU +
       " RGSM=" + RGSM +
       conditionalParameter(RGCN != null && !RGCN.isEmpty, " RGCN=" + RGCN) +
       conditionalParameter(RGDS != null && !RGDS.isEmpty, " RGDS=" + RGDS)
}