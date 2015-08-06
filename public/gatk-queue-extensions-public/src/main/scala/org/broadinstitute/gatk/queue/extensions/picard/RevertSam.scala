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
class RevertSam extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "RevertSam"
  javaMainClass = "picard.sam.RevertSam"

  @Input(shortName = "input", fullName = "input_bam_files", required = true, doc = "The input SAM or BAM files to revert.")
  var input: Seq[File] = Nil

  @Output(shortName = "output", fullName = "output_bam_file", required = true, doc = "The reverted BAM or SAM output file.")
  var output: File = _

  @Output(shortName = "out_index", fullName = "output_bam_index_file", required = false, doc = "The output bam index")
  var outputIndex: File = _

  @Argument(shortName = "roq", fullName = "restore_original_qualities", required = false, doc = "True to restore original qualities from the OQ field to the QUAL field if available.")
  var restoreOriginalQualities: Boolean = true

  @Argument(shortName = "rdi", fullName = "remove_duplicate_information", required = false, doc = "Remove duplicate read flags from all reads. Note that if this is true and REMOVE_ALIGNMENT_INFORMATION==false, the output may have the unusual but sometimes desirable trait of having unmapped reads that are marked as duplicates.")
  var removeDuplicateInformation: Boolean = true

  @Argument(shortName = "rai", fullName = "remove_alignment_information", required = false, doc = "Remove all alignment information from the file.")
  var removeAlignmentInformation: Boolean = true

  @Argument(shortName = "atc", fullName = "attributes_to_clear", required = false, doc = "When removing alignment information, the set of optional tags to remove.")
  var attributesToClear: Seq[String] = Nil

  @Argument(shortName = "sa", fullName = "sample_alias", required = false, doc = "The sample alias to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias.")
  var sampleAlias: String = null

  @Argument(shortName = "ln", fullName = "library_name", required = false, doc = "The library name to use in the reverted output file. This will override the existing sample alias in the file and is used only if all the read groups in the input file have the same sample alias.")
  var libraryName: String = null

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (outputIndex == null && output != null)
      outputIndex = new File(output.getName.stripSuffix(".bam") + ".bai")
  }


  override def inputBams = input
  override def outputBam = output
  this.createIndex = Some(true)
  override def commandLine = super.commandLine +
                             conditional(!restoreOriginalQualities, "RESTORE_ORIGINAL_QUALITIES=false") +
                             conditional(!removeDuplicateInformation, "REMOVE_DUPLICATE_INFORMATION=false") +
                             conditional(!removeAlignmentInformation, "REMOVE_ALIGNMENT_INFORMATION=false") +
                             repeat("ATTRIBUTE_TO_CLEAR=", attributesToClear, spaceSeparated=false) +   // repeat() returns "" for null/empty list
                             conditional(sampleAlias != null, "SAMPLE_ALIAS=" + sampleAlias) +
                             conditional(libraryName != null, "LIBRARY_NAME=" + libraryName)
}