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

package org.broadinstitute.gatk.queue.extensions.samtools

import java.io.File
import org.broadinstitute.gatk.utils.commandline.{Argument, Output, Input}

/**
 * Merges BAM files using samtools.
 */
class SamtoolsMergeFunction extends SamtoolsCommandLineFunction {
  analysisName = "samtools merge"

  @Input(doc="BAM file input")
  var inputBams: Seq[File] = Nil

  @Output(doc="BAM file output")
  var outputBam: File = _

  @Argument(doc="region", required=false)
  var region: String = _

  @Input(doc="BAM file input indexes")
  var inputBamIndexes: Seq[File] = Nil

  override def freezeFieldValues() {
    super.freezeFieldValues()
    inputBamIndexes ++= inputBams
      .filter(orig => orig != null && orig.getName.endsWith(".bam"))
      .flatMap(orig => Array(
        new File(orig.getPath + ".bai"),
        new File(orig.getPath.stripSuffix(".bam") + ".bai")
      ))
  }

  def commandLine = required(samtools) +
                    required("merge") +
                    optional("-R", region) +
                    required(outputBam) +
                    repeat(inputBams)
}
