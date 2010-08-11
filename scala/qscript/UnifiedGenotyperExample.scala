import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class UnifiedGenotyperExample extends QScript {
  qscript =>

  @Input(doc="gatk jar file")
  var gatkJar: File = _

  @Input(doc="bam files", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="interval list", shortName="L")
  var intervals: File = _

  @Input(doc="referenceFile", shortName="R")
  var referenceFile: File = _

  @Argument(doc="filter names", shortName="filter")
  var filterNames: List[String] = Nil

  @Argument(doc="filter expressions", shortName="filterExpression")
  var filterExpressions: List[String] = Nil

  @Argument(doc="job queue", shortName="queue", required=false)
  var jobQueue = "broad"

  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.jobQueue = qscript.jobQueue
    this.jarFile = qscript.gatkJar
    this.intervals = qscript.intervals
    this.reference_sequence = qscript.referenceFile
  }

  def script = {
    for (bam <- bamFiles) {
      val ug = new UnifiedGenotyper with UnifiedGenotyperArguments
      val vf = new VariantFiltration with UnifiedGenotyperArguments
      val ve = new VariantEval with UnifiedGenotyperArguments

      // Make sure the Sting/shell folder is in your path to use mergeText.sh and splitIntervals.sh.
      ug.scatterCount = 3
      ug.cleanupTempDirectories = true
      ug.input_file :+= bam
      ug.out = swapExt(bam, "bam", "unfiltered.vcf")

      vf.rodBind :+= RodBind("vcf", "VCF", ug.out)
      vf.out = swapExt(bam, "bam", "filtered.vcf")

      ve.rodBind :+= RodBind("vcf", "VCF", vf.out)
      ve.out = swapExt(bam, "bam", "eval")

      add(ug, vf, ve)
    }
  }
}
