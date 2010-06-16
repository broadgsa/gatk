package org.broadinstitute.sting.queue.function.gatk

import java.io.File
import org.broadinstitute.sting.queue.util.{Input, Optional}
import org.broadinstitute.sting.queue.function.{MemoryLimitedFunction, IntervalFunction, CommandLineFunction}

trait GatkFunction extends CommandLineFunction with MemoryLimitedFunction with IntervalFunction {
  @Input
  @Optional
  var javaTmpDir: String = _

  @Input
  var gatkJar: String = _

  @Input
  var referenceFile: String = _

  @Input
  @Optional
  var bamFiles: List[File] = Nil

  @Input
  @Optional
  var dbsnp: File = _

  @Input
  @Optional
  var intervals: Intervals = new Intervals("all")

  protected def gatkCommandLine(walker: String) =
    "java%s%s -jar %s -T %s -R %s%s%s "
    .format(optional(" -Xmx", memoryLimit, "g"), optional(" -Djava.io.tmpdir=", javaTmpDir),
      gatkJar, walker, referenceFile, repeat(" -I ", bamFiles), optional(" -D ", dbsnp), optional(" -L ", intervals))
}
