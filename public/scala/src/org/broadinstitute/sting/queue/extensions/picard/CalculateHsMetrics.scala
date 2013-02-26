/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import java.io.File
import net.sf.picard.analysis.MetricAccumulationLevel

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 10/9/12
 * Time: 5:59 PM
 * To change this template use File | Settings | File Templates.
 */
class CalculateHsMetrics extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardMetricsFunction {
  analysisName = "CalculateHsMetrics"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

  @Output(doc="The output file to write statistics to", shortName = "output", fullName = "output_file", required = true)
  var output: File = _

  @Argument(doc="Interval list with targets", shortName = "targets", fullName = "target_list", required = true)
  var targets: File = _

  @Argument(doc="Interval list with baits", shortName = "baits", fullName = "bait_list", required = true)
  var baits: File = _

  @Argument(doc="Reference file", shortName = "reference", fullName = "reference", required = true)
  var reference: File = _

  val level = MetricAccumulationLevel.SAMPLE

  override def inputBams = input
  override def outputFile = output
  override def commandLine = super.commandLine +
    required("BAIT_INTERVALS=" + baits) +
    required("TARGET_INTERVALS=" + targets) +
    required("REFERENCE_SEQUENCE=" + reference) +
    optional("METRIC_ACCUMULATION_LEVEL="+level)

}
