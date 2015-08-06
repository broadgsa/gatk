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
class AddOrReplaceReadGroups extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "AddOrReplaceReadGroups"
  javaMainClass = "picard.sam.AddOrReplaceReadGroups"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

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
                             required("RGID=" + RGID) +
                             required("RGLB=" + RGLB) +
                             required("RGPL=" + RGPL) +
                             required("RGPU=" + RGPU) +
                             required("RGSM=" + RGSM) +
                             conditional(RGCN != null && !RGCN.isEmpty, "RGCN=" + RGCN) +
                             conditional(RGDS != null && !RGDS.isEmpty, "RGDS=" + RGDS)
}