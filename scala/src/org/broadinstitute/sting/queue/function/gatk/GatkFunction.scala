package org.broadinstitute.sting.queue.function.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.IntervalFunction
import org.broadinstitute.sting.queue.function.scattergather.{Scatter, ScatterGatherableFunction, IntervalScatterFunction}
import org.broadinstitute.sting.commandline.Input

trait GatkFunction extends ScatterGatherableFunction with IntervalFunction {
  @Input(doc="Temporary directory to write any files", required=false)
  var javaTmpDir: String = _

  @Input(doc="GATK jar")
  var gatkJar: String = _

  @Input(doc="Reference fasta")
  var referenceFile: File = _

  @Input(doc="Bam files", required=false)
  var bamFiles: List[File] = Nil

  @Input(doc="Intervals", required=false)
  @Scatter(classOf[IntervalScatterFunction])
  var intervals: File = _

  @Input(doc="DBSNP", required=false)
  var dbsnp: File = _

  protected def gatkCommandLine(walker: String) =
    "java%s%s -jar %s -T %s -R %s%s%s%s "
    .format(optional(" -Xmx", memoryLimit, "g"), optional(" -Djava.io.tmpdir=", javaTmpDir),
      gatkJar, walker, referenceFile, repeat(" -I ", bamFiles), optional(" -D ", dbsnp), optional(" -L ", intervals))
}
