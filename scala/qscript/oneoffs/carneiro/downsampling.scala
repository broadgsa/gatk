package oneoffs.carneiro

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import scala.io.Source._
/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 3/17/11
 * Time: 11:29 AM
 * To change this template use File | Settings | File Templates.
 */


class downsampling extends QScript {

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="target intervals", shortName="t", required=true)
  var targetIntervals: File = _

  @Input(doc="bootstrap number", shortName="b", required=false)
  var bootstrap: Int = 1

  @Input(doc="downsampling step", shortName="ds", required=true)
  var downsamplingStep: Double = _

  @Input(doc="downsampling floor", shortName="df", required=false)
  var downsamplingFloor: Double = 0.0

  @Input(doc="downsampling ceiling", shortName="dc", required=false)
  var downsamplingCeiling: Double = 1.0

  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbSNP file", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf")

  @Input(doc="project name", shortName="p", required=false)
  var base: String = "prj"

  def countLines(file: File):Int = {
    var count: Int = 0
    for (l <- fromFile(file).getLines) {
      count = count + 1
    }
    return count
  }

  val queueLogDir: String = ".qlog/"
  val outFile: String = "cov.out"
  val fullCoverageVCF = new File("/humgen/gsa-hpprojects/dev/carneiro/downsampling/analysis/fullcov/fullcov.F1.filtered.vcf")

  def script = {
    val nIntervals = math.min(200, countLines(targetIntervals))

    var f: Double = downsamplingCeiling
    var i: Int = 1
    while (f>=downsamplingFloor) {
      var b: Int = bootstrap
      while(b > 0) {
        val file = swapExt(outFile, ".out", ".F" + i + "." + b + ".out")
        add(cov(f, file))
        b = b - 1
      }
      val snp_out = new File(base + ".F" + i + ".raw.vcf")
      val filter_out = new File(base + ".F" + i + ".filtered.vcf")
      val eval_out = new File(base + ".F" + i + ".eval")

      add( snps(f, snp_out, nIntervals),
           filter(snp_out, filter_out),
           eval(filter_out, eval_out))

      f = f - downsamplingStep
      i = i + 1
    }
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals :+= targetIntervals
    this.jarFile = GATKjar
    this.reference_sequence = reference
    this.memoryLimit = 4
  }

  case class cov (fraction: Double, outFile: File) extends Percent20xCoverage with CommandLineGATKArgs {
    this.input_file :+= input
    this.out = outFile
    this.ndrs = true
    this.downsample_to_fraction = fraction
    this.jobName = queueLogDir + outFile + ".cov"
  }

  case class snps (fraction: Double, outFile: File, nIntervals: Int) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.memoryLimit = 6
    this.downsample_to_coverage = 600
    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.input_file :+= input
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.downsample_to_fraction = fraction
    this.scatterCount = nIntervals
    this.out = outFile
    this.analysisName = outFile + "_snps"
    this.jobName = queueLogDir + outFile
  }

  case class filter (inFile: File, outFile: File) extends VariantFiltration with CommandLineGATKArgs {
    this.filterName ++= List("SNPSBFilter","SNPQDFilter","SNPHRunFilter")
    this.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    this.clusterWindowSize = 10
    this.clusterSize = 3
    this.variantVCF = inFile
    this.out = outFile
    this.analysisName = outFile + "_filter"
    this.jobName = queueLogDir + outFile
  }

  case class eval (inFile: File, outFile: File) extends VariantEval with CommandLineGATKArgs {
    this.noST = true
    this.noEV = true
    this.evalModule ++= List("TiTvVariantEvaluator", "CountVariants", "ValidationReport")
    this.stratificationModule ++= List("EvalRod", "CompRod", "Novelty")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.rodBind :+= RodBind("eval", "VCF", inFile)
    this.rodBind :+= RodBind("comp", "VCF", fullCoverageVCF)
    this.out = outFile
    this.analysisName = outFile + "_VariantEval"
    this.jobName = queueLogDir + outFile
  }
}
