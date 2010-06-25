import org.broadinstitute.sting.queue.QScript._

setArgs(args)

for (bam <- inputs("bam")) {
  val ug = new UnifiedGenotyper
  val vf = new VariantFiltration
  val ve = new GatkFunction {
    @Input(doc="vcf") var vcfFile: File = _
    @Output(doc="eval") var evalFile: File = _
    def commandLine = gatkCommandLine("VariantEval") + "-B eval,VCF,%s -o %s".format(vcfFile, evalFile)
  }

  // Make sure the Sting/shell folder is in your path to use mergeText.sh and splitIntervals.sh.
  ug.scatterCount = 3
  ug.bamFiles :+= bam
  ug.vcfFile = swapExt(bam, "bam", "unfiltered.vcf")

  vf.vcfInput = ug.vcfFile
  vf.vcfOutput = swapExt(bam, "bam", "filtered.vcf")

  ve.vcfFile = vf.vcfOutput
  ve.evalFile = swapExt(bam, "bam", "eval")

  add(ug, vf, ve)
}

setParams
run


class UnifiedGenotyper extends GatkFunction {
  @Output(doc="vcf")
  @Gather(classOf[SimpleTextGatherFunction])
  var vcfFile: File = _
  def commandLine = gatkCommandLine("UnifiedGenotyper") + "-varout %s".format(vcfFile)
}

class VariantFiltration extends GatkFunction {
  @Input(doc="input vcf")
  var vcfInput: File = _

  @Input(doc="filter names")
  var filterNames: List[String] = Nil

  @Input(doc="filter expressions")
  var filterExpressions: List[String] = Nil

  @Output(doc="output vcf")
  var vcfOutput: File = _

  def commandLine = gatkCommandLine("VariantFiltration") + "%s%s -B variant,VCF,%s -o %s"
      .format(repeat(" -filterName ", filterNames), repeat(" -filterExpression ", filterExpressions), vcfInput, vcfOutput)
}
