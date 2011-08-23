package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline._

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class RevertSam extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "RevertSam"
  javaMainClass = "net.sf.picard.sam.RevertSam"

  @Input(shortName = "input", fullName = "input_bam_files", required = true, doc = "The input SAM or BAM files to revert.")
  var input: List[File] = Nil

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
  var attributesToClear: List[String] = Nil

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
       conditionalParameter(!restoreOriginalQualities, " RESTORE_ORIGINAL_QUALITIES=false") +
       conditionalParameter(!removeDuplicateInformation, " REMOVE_DUPLICATE_INFORMATION=false") +
       conditionalParameter(!removeAlignmentInformation, " REMOVE_ALIGNMENT_INFORMATION=false") +
       conditionalParameter(!attributesToClear.isEmpty, repeat(" ATTRIBUTE_TO_CLEAR=", attributesToClear)) +
       conditionalParameter(sampleAlias != null, " SAMPLE_ALIAS=" + sampleAlias) +
       conditionalParameter(libraryName != null, " LIBRARY_NAME=" + libraryName)
}