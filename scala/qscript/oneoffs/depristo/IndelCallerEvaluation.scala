package oneoffs.depristo

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction

class IndelCallerEvaluation extends QScript {
  val BUNDLE = "/humgen/gsa-hpprojects/GATK/bundle/current"

  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File(BUNDLE + "/b37/human_g1k_v37.fasta")

  @Argument(shortName = "bam", doc="BAM", required=true)
  val bams: List[File] = null;

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val myIntervals: String = null;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 250;

  val dbSNP: File = new File(BUNDLE + "/b37/dbsnp_132.b37.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4

    if ( intervals != null )
      this.intervalsString = List(myIntervals);
  }

  trait CoFoJa extends JavaCommandLineFunction {
    override def javaOpts = super.javaOpts //  + " -javaagent:lib/cofoja.jar"
  }

  def processOne(bam: File, gsaProduction: Boolean): File = {
    val rawVCF = new Call(bam, gsaProduction)
    add(rawVCF)

    val filterIndels = new FilterIndels(rawVCF.out)
    add(filterIndels)

    // create a variant eval for us
    add(new Eval(filterIndels.out))
    return filterIndels.out
  }

  def script = {
    for ( gsaProduction <- List(true, false)) {
      val vcfs = bams.map(processOne(_, gsaProduction))

      val combineCalls = new CombineVariants with UNIVERSAL_GATK_ARGS
      for ( vcf <- vcfs )
        combineCalls.rodBind :+= RodBind(vcf.getName, "VCF", vcf)

      combineCalls.filteredrecordsmergetype = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
      combineCalls.out = "combined" + productionString(gsaProduction) + ".vcf"
      add(combineCalls)

      add(new ToTable(combineCalls.out))
    }
  }

  class FilterIndels(@Input vcf: File) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.variantVCF = vcf
    this.filterName = List("Indel_QUAL", "Indel_SB", "Indel_QD")
    this.filterExpression = List("\"QUAL<30.0\"", "\"SB>-1.0\"", "\"QD<2.0\"")
    this.out = swapExt(vcf,".vcf",".filtered.vcf")
  }

  class ToTable(@Input vcf: File) extends VariantsToTable with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("variant", "VCF", vcf)
    this.fields = List("FILTER", "set")
    this.out = swapExt(vcf,".vcf",".table")
    this.raw = true
  }

  class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.evalModule = List("CountVariants", "IndelStatistics", "CompOverlap")
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "JexlExpression")
    this.out = swapExt(vcf,".vcf",".eval")
  }

  def productionString(gsaProduction: Boolean): String = {
    return if ( gsaProduction ) ".prod" else ".expt"
  }

  class Call(@Input(doc="foo") bam: File, gsaProduction: Boolean) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = swapExt(bam,".bam", productionString(gsaProduction) + ".indels.vcf")
    this.input_file = List(bam)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.o = outVCF

    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.GSA_PRODUCTION_ONLY = gsaProduction

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
  }
}

