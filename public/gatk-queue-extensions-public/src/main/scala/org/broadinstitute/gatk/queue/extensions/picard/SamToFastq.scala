/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.extensions.picard

import org.broadinstitute.gatk.utils.commandline._

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class SamToFastq extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "SamToFastq"
  javaMainClass = "picard.sam.SamToFastq"

  @Input(shortName = "input", fullName = "input_bam_files", required = true, doc = "Input SAM/BAM file to extract reads from.")
  var input: Seq[File] = Nil

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

  @Argument(shortName = "il", fullName = "interleave", required = false, doc = "Will generate an interleaved fastq if paired, each line will have /1 or /2 to describe which end it came from")
  var interleave: Boolean = false

  override def inputBams = input
  override def outputBam = null
  this.sortOrder = null

  override def commandLine = super.commandLine +
                             required("FASTQ=" + fastq) +
                             optional("SECOND_END_FASTQ=", secondEndFastQ, spaceSeparated=false) +
                             conditional(outputPerReadGroup, "OUTPUT_PER_RG=" + outputPerReadGroup) +
                             optional("OUTPUT_DIR=", outputDir, spaceSeparated=false) +
                             conditional(!reReverse, "RE_REVERSE=" + reReverse) +
                             conditional(includeNonPFReads, "INCLUDE_NON_PF_READS=" + includeNonPFReads) +
                             optional("CLIPPING_ATTRIBUTE=", clippingAttribute, spaceSeparated=false) +
                             optional("CLIPPING_ACTION=", clippingAction, spaceSeparated=false) +
                             conditional(readOneTrim >= 0, "READ1_TRIM=" + readOneTrim) +
                             conditional(readOneMaxBasesToWrite >= 0, "READ1_MAX_BASES_TO_WRITE=" + readOneMaxBasesToWrite) +
                             conditional(readTwoTrim >= 0, "READ2_TRIM=" + readTwoTrim) +
                             conditional(readTwoMaxBasesToWrite >= 0, "READ2_MAX_BASES_TO_WRITE=" + readTwoMaxBasesToWrite) +
                             conditional(includeNonPrimaryAlignments, "INCLUDE_NON_PRIMARY_ALIGNMENTS=" + includeNonPrimaryAlignments) +
                             conditional(interleave, "INTERLEAVE=" + interleave)

}