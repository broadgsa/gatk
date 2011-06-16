package oneoffs.depristo

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction

class ReducedBAMEvaluation extends QScript {
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "bam", doc="BAM", required=true)
  val bam: File = null;

  @Argument(shortName = "reduceIntervals", doc="intervals", required=false)
  val REDUCE_INTERVAL: String = null;

  @Argument(shortName = "callingIntervals", doc="intervals", required=false)
  val CALLING_INTERVAL: String = null;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 250;

  @Argument(shortName = "minimalVCF", doc="", required=false)
  val minimalVCF: Boolean = false;

  val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_132.b37.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }

  trait CoFoJa extends JavaCommandLineFunction {
    override def javaOpts = super.javaOpts //  + " -javaagent:lib/cofoja.jar"
  }

  def script = {
    val reducedBAM = new ReduceBAM(bam)
    add(reducedBAM)
    val reducedBAMVCF = callAndEvaluateBAM(reducedBAM.out)

    val slicedBAM = new SliceBAM(bam)
    add(slicedBAM)
    val fullBAMVCF = callAndEvaluateBAM(slicedBAM.out)

    val combineCalls = new CombineVariants with UNIVERSAL_GATK_ARGS
    combineCalls.rodBind :+= RodBind("fullBAM", "VCF", fullBAMVCF)
    combineCalls.rodBind :+= RodBind("reducedBAM", "VCF", reducedBAMVCF)
    combineCalls.rod_priority_list = "reducedBAM,fullBAM"
    combineCalls.filteredrecordsmergetype = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineCalls.out = swapExt(reducedBAM.out,".bam",".filtered.combined.vcf")
    add(combineCalls)
    val eval = new Eval(combineCalls.out) // evaluate the combined VCF
    eval.select = List("'set==\"Intersection\"'", "'set==\"fullBAM\"'", "'set==\"reducedBAM\"'", "'set==\"filterInreducedBAM-fullBAM\"'", "'set==\"reducedBAM-filterInfullBAM\"'")
    eval.selectName = List("Intersection", "fullBAM", "reducedBAM", "filterInreducedBAM-fullBAM", "reducedBAM-filterInfullBAM")
    add(eval)
  }

  def callAndEvaluateBAM(bam: File): File = {
    val rawVCF = new Call(bam)
    add(rawVCF)

    val filterSNPs = new VariantFiltration with UNIVERSAL_GATK_ARGS
    filterSNPs.variantVCF = rawVCF.out
    filterSNPs.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    filterSNPs.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    filterSNPs.clusterWindowSize = 10
    filterSNPs.clusterSize = 3
    filterSNPs.out = swapExt(rawVCF.out,".vcf",".filtered.vcf")
    add(filterSNPs)

    // create a variant eval for us
    add(new Eval(filterSNPs.out))

    // for convenient diffing
    add(new DiffableTable(rawVCF.out))
    add(new DiffableTable(filterSNPs.out))

    return filterSNPs.out
  }

  class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants")
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "JexlExpression")
    this.out = swapExt(vcf,".vcf",".eval")
    if ( CALLING_INTERVAL != null )
      this.intervalsString = List(CALLING_INTERVAL);
  }


  class ReduceBAM(bam: File) extends ReduceReads with UNIVERSAL_GATK_ARGS with CoFoJa {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = swapExt(bam,".bam",".reduced.bam")
    this.CS = 20
    this.mravs = 50
    this.mbrc = 10000

    if ( REDUCE_INTERVAL != null )
      this.intervalsString = List(REDUCE_INTERVAL);
  }

  class SliceBAM(bam: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = swapExt(bam,".bam",".printreads.bam")
    if ( REDUCE_INTERVAL != null )
      this.intervalsString = List(REDUCE_INTERVAL);
  }

  class Call(@Input(doc="foo") bam: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = swapExt(bam,".bam",".vcf")
    this.input_file = List(bam)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.o = outVCF

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    if ( minimalVCF )
      this.group = List("none")

    if ( CALLING_INTERVAL != null ) {
      this.intervalsString = List(CALLING_INTERVAL)
    }
  }

  class DiffableTable(@Input vcf: File) extends CommandLineFunction {
    @Output var out: File = swapExt(vcf,".vcf",".table")
    def commandLine = "cut -f 1,2,4,5,7 %s > %s".format(vcf, out)
  }
}

