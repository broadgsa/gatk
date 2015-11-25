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

import org.broadinstitute.gatk.utils.commandline.{Argument, Output, Input}
import java.io.File

class CollectWgsMetrics extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardMetricsFunction {
  analysisName = "CollectWgsMetrics"
  javaMainClass = "picard.analysis.CollectWgsMetrics"

  @Input(doc = "The input SAM or BAM files to analyze", shortName = "i", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

  @Output(doc = "The output file to write statistics to", shortName = "o", fullName = "output_file", required = true)
  var output: File = _

  @Argument(doc = "Reference file", shortName = "r", fullName = "reference", required = true)
  var reference: File = _

  @Argument(doc = "Minimum mapping quality for a read to contribute coverage.", shortName = "mq", fullName = "minimum_mapping_quality", required = false)
  var mq: Integer = _

  @Argument(doc = "Minimum base quality for a base to contribute coverage.", shortName = "q", fullName = "minimum_base_quality", required = false)
  var q: Integer = _

  @Argument(doc = "Treat bases with coverage exceeding this value as if they had coverage at this value.", shortName = "cap", fullName = "coverage_cap", required = false)
  var cap: Integer = _

  @Argument(doc = "For debugging purposes, stop after processing this many genomic bases.", fullName = "stop_after", required = false)
  var stopAfter: Long = _

  @Argument(doc = "Determines whether to include the base quality histogram in the metrics file.", fullName = "include_bq_histogram", required = false)
  var includeBQHistogram: Boolean = _

  override def inputBams = input

  override def outputFile = output

  override def commandLine = super.commandLine +
    required("REFERENCE_SEQUENCE=" + reference) +
    optional("MQ=", mq, spaceSeparated = false) +
    optional("Q=", q, spaceSeparated = false) +
    optional("CAP=", cap, spaceSeparated = false) +
    optional("STOP_AFTER=", stopAfter, spaceSeparated = false) +
    optional("INCLUDE_BQ_HISTOGRAM=", includeBQHistogram, spaceSeparated = false)
}