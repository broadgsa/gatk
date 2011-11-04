package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class SamToFastq extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "SamToFastq"
  javaMainClass = "net.sf.picard.sam.SamToFastq"

  @Input(shortName = "input", fullName = "input_bam_files", required = true, doc = "Input SAM/BAM file to extract reads from.")
  var input: List[File] = Nil

  @Output(shortName = "fastq", fullName = "output_fastq_file", required = true, doc = "Output fastq file (single-end fastq or, if paired, first end of the pair fastq).")
  var fastq: File = _

  @Output(shortName = "se", fullName = "second_end_fastq", required = false, doc = "Output fastq file (if paired, second end of the pair fastq).")
  var secondEndFastQ: File = _

  @Argument(shortName = "opg", fullName = "output_per_readgroup", required = false, doc = "Output a fastq file per read group (two fastq files per read group if the group is paired).")
  var outputPerReadGroup: Boolean = false

  @Argument(shortName = "od", fullName = "output_dir", required = false, doc = "Directory in which to output the fastq file(s). Used only when OUTPUT_PER_RG is true.")
  var outputDir: File = _

  @Argument(shortName = "rr", fullName = "re_reverse", required = false, doc = "Re-reverse bases and qualities of reads with negative strand flag set before writing them to fastq.")
  var reReverse: Boolean = true

  @Argument(shortName = "nonpf", fullName = "include_non_pf_reads", required = false, doc = "Include non-PF reads from the SAM file into the output FASTQ files.")
  var includeNonPFReads: Boolean = false

  @Argument(shortName = "cat", fullName = "clipping_attribute", required = false, doc = "The attribute that stores the position at which the SAM record should be clipped.")
  var clippingAttribute: String = null

  @Argument(shortName = "cac", fullName = "clipping_action", required = false, doc = "The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position; 'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region.")
  var clippingAction: String = null

  @Argument(shortName = "r1t", fullName = "read_one_trim", required = false, doc = "The number of bases to trim from the beginning of read 1.")
  var readOneTrim: Int = -1

  @Argument(shortName = "r1mbtw", fullName = "read_one_max_bases_to_write", required = false, doc = "The maximum number of bases to write from read 1 after trimming. If there are fewer than this many bases left after trimming, all will be written. If this value is null then all bases left after trimming will be written.")
  var readOneMaxBasesToWrite: Int = -1

  @Argument(shortName = "r2t", fullName = "read_two_trim", required = false, doc = "The number of bases to trim from the beginning of read 2.")
  var readTwoTrim: Int = -1

  @Argument(shortName = "r2mbtw", fullName = "read_two_max_bases_to_write", required = false, doc = "The maximum number of bases to write from read 2 after trimming. If there are fewer than this many bases left after trimming, all will be written. If this value is null then all bases left after trimming will be written.")
  var readTwoMaxBasesToWrite: Int = -1

  @Argument(shortName = "inpa", fullName = "include_non_primary_alignments", required = false, doc = "If true, include non-primary alignments in the output. Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.")
  var includeNonPrimaryAlignments: Boolean = false

  override def inputBams = input
  override def outputBam = null
  this.sortOrder = null

  override def commandLine = super.commandLine +
       " FASTQ=" + fastq +
       optional(" SECOND_END_FASTQ=", secondEndFastQ) +
       conditionalParameter(outputPerReadGroup, optional(" OUTPUT_PER_RG=", outputPerReadGroup)) +
       optional(" OUTPUT_DIR=", outputDir) +
       conditionalParameter(!reReverse, optional(" RE_REVERSE=", reReverse)) +
       conditionalParameter(includeNonPFReads, optional(" INCLUDE_NON_PF_READS=", includeNonPFReads)) +
       optional(" CLIPPING_ATTRIBUTE=", clippingAttribute) +
       optional(" CLIPPING_ACTION=", clippingAction) +
       conditionalParameter (readOneTrim >= 0, optional(" READ1_TRIM=", readOneTrim)) +
       conditionalParameter (readOneMaxBasesToWrite >= 0, optional(" READ1_MAX_BASES_TO_WRITE=", readOneMaxBasesToWrite)) +
       conditionalParameter (readTwoTrim >= 0, optional(" READ2_TRIM=", readTwoTrim)) +
       conditionalParameter (readTwoMaxBasesToWrite >=0, optional(" READ2_MAX_BASES_TO_WRITE=", readTwoMaxBasesToWrite)) +
       conditionalParameter (includeNonPrimaryAlignments, optional(" INCLUDE_NON_PRIMARY_ALIGNMENTS=", includeNonPrimaryAlignments))
}