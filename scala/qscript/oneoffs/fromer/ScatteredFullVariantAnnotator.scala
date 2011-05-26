package oneoffs.fromer

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.interval.IntervalSetRule

class ScatteredFullVariantAnnotator extends QScript {
  qscript =>

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(shortName = "L", doc = "Intervals", required = false)
  var intervals: String = null

  @Input(doc = "level of parallelism. By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCount = 0

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "variant calls to annotate", fullName = "variantVCF", shortName = "C", required = true)
  var variantVCF: File = _

  @Output(doc = "annotated file to output", shortName = "o", required = true)
  var outputAnnotated: File = _

  @Output(doc = "Memory limit", fullName = "memoryLimit", shortName = "m", required = false)
  var memoryLimit = 3

  def script = {
    add(new ScatteredFullVariantAnnotator())
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    if (qscript.intervals != null) {
      this.intervalsString = List(qscript.intervals)
    }
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.input_file = List(qscript.bams)

    this.memoryLimit = qscript.memoryLimit
    this.logging_level = "INFO"

    this.rodToIntervalTrackName = "variant"
    this.BTI_merge_rule = IntervalSetRule.INTERSECTION
  }

  class ScatteredFullVariantAnnotator() extends org.broadinstitute.sting.queue.extensions.gatk.VariantAnnotator with CommandLineGATKArgs {
    this.scatterCount = qscript.scatterCount
    this.variantVCF = qscript.variantVCF
    this.useAllAnnotations = true
    this.out = qscript.outputAnnotated
  }
}