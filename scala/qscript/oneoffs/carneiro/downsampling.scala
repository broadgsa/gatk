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

  @Input(doc="HapMap file", shortName="H", required=false)
  var hapmap: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf")

  @Input(doc="Omni file", shortName="O", required=false)
  var omni: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf")

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
  val trancheTarget = "99.0"

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

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  case class VQSR(inFile: File, tranchesFiles: File, outFile: File) extends VariantRecalibrator with CommandLineGATKArgs {
    this.rodBind :+= RodBind("input", "VCF", inFile)
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni, "known=false,training=true,truth=true,prior=12.0")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP, "known=true,training=false,truth=false,prior=10.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "HRun")
    this.tranches_file = tranchesFile
    this.recal_file = outFile
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = t.name + "_VQSR"
    this.jobName =  queueLogDir + outFile
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (inFile: File, tranchesFiles: File, outFile: File) extends ApplyRecalibration with CommandLineGATKArgs {
    this.rodBind :+= RodBind("input", "VCF", inFile)
    this.tranches_file = tranchesFile
    this.recal_file = inFile
    this.ts_filter_level = trancheTarget
    this.out = outFile
    this.analysisName = outFile + "_AVQSR"
    this.jobName =  queueLogDir + outFile
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
