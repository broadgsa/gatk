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
