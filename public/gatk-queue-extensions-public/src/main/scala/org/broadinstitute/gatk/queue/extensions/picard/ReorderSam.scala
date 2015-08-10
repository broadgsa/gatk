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
class ReorderSam extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "ReorderSam"
  javaMainClass = "picard.sam.ReorderSam"

  @Input(doc="Input file (bam or sam) to extract reads from.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

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
                             required("REFERENCE=" + sortReference) +
                             optional("ALLOW_INCOMPLETE_DICT_CONCORDANCE=", ALLOW_INCOMPLETE_DICT_CONCORDANCE, spaceSeparated=false)
                             optional("ALLOW_CONTIG_LENGTH_DISCORDANCE=", ALLOW_CONTIG_LENGTH_DISCORDANCE, spaceSeparated=false)
}