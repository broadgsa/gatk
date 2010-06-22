package org.broadinstitute.sting.queue.function.gatk

import java.io.File
import org.broadinstitute.sting.queue.function.IntervalFunction
import org.broadinstitute.sting.queue.util.{Scatter, Internal, Input, Optional}
import org.broadinstitute.sting.queue.function.scattergather.{ScatterGatherableFunction, IntervalScatterFunction}

trait GatkFunction extends ScatterGatherableFunction with IntervalFunction {
  @Internal
  @Optional
  var javaTmpDir: String = _

  @Input
  var gatkJar: String = _

  @Input
  var referenceFile: File = _

  @Input
  @Optional
  var bamFiles: List[File] = Nil

  @Input
  @Optional
  @Scatter(classOf[IntervalScatterFunction])
  var intervals: File = _

  @Input
  @Optional
  var dbsnp: File = _

  protected def gatkCommandLine(walker: String) =
    "java%s%s -jar %s -T %s -R %s%s%s%s "
    .format(optional(" -Xmx", memoryLimit, "g"), optional(" -Djava.io.tmpdir=", javaTmpDir),
      gatkJar, walker, referenceFile, repeat(" -I ", bamFiles), optional(" -D ", dbsnp), optional(" -L ", intervals))
}
