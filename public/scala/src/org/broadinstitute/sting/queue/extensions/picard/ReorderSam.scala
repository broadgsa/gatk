package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File
/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class ReorderSam extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "ReorderSam"
  javaMainClass = "net.sf.picard.sam.ReorderSam"

  @Input(doc="Input file (bam or sam) to extract reads from.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: List[File] = Nil

  @Output(doc="Output file (bam or sam) to write extracted reads to.", shortName = "output", fullName = "output_bam_file", required = true)
  var output: File = _

  @Output(doc="The output bam index", shortName = "out_index", fullName = "output_bam_index_file", required = false)
  var outputIndex: File = _

  @Argument(doc="Reference sequence to reorder reads to match.", shortName = "ref", fullName = "sort_reference", required = true)
  var sortReference: File = _

  @Argument(doc="If true, then allows only a partial overlap of the BAM contigs with the new reference sequence contigs. By default, this tool requires a corresponding contig in the new reference for each read contig.", shortName = "aic", fullName = "allow_incomplete_concordance", required = false)
  var ALLOW_INCOMPLETE_DICT_CONCORDANCE: Boolean = _

  @Argument(doc="If true, then permits mapping from a read contig to a new reference contig with the same name but a different length. Highly dangerous, only use if you know what you are doing.", shortName = "acld", fullName = "allow_contig_length_discordance", required = false)
  var ALLOW_CONTIG_LENGTH_DISCORDANCE: Boolean = _

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }

  override def inputBams = input
  override def outputBam = output
  this.createIndex = Some(true)
  this.sortOrder = null
  override def commandLine = super.commandLine +
    " REFERENCE=" + sortReference +
    optional(" ALLOW_INCOMPLETE_DICT_CONCORDANCE=", ALLOW_INCOMPLETE_DICT_CONCORDANCE)
    optional(" ALLOW_CONTIG_LENGTH_DISCORDANCE=", ALLOW_CONTIG_LENGTH_DISCORDANCE)
}