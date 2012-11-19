package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import java.io.File

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 10/10/12
 * Time: 10:37 AM
 * To change this template use File | Settings | File Templates.
 */
class CollectGcBiasMetrics extends org.broadinstitute.sting.queue.function.JavaCommandLineFunction with PicardMetricsFunction {
  analysisName = "CollectGcBiasMetrics"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

  @Output(doc="The output file to write statistics to", shortName = "output", fullName = "output_file", required = true)
  var output: File = _

  @Argument(doc="Reference file", shortName = "reference", fullName = "reference", required = true)
  var reference: File = _

  override def inputBams = input
  override def outputFile = output
  override def commandLine = super.commandLine +
    required("SUMMARY_OUTPUT=" + output) +
    required("CHART_OUTPUT=" + output+".pdf") +
    required("REFERENCE_SEQUENCE=" + reference) +
    required("ASSUME_SORTED=true")
}
