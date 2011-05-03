package oneoffs.carneiro

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class justRecalibrateEach extends QScript {

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to AnalyzeCovariates.jar", fullName="path_to_ac_jar", shortName="ac", required=false)
  var ACJar: File = new File("/humgen/gsa-scr1/carneiro/stable/dist/AnalyzeCovariates.jar")

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=false)
  var R: String = new File("/humgen/gsa-scr1/carneiro/stable/R")

  @Input(doc="bad regions interval", shortName="bad", required=false)
  var badInterval: File = new File("/humgen/gsa-hpprojects/dev/carneiro/goodbad/data/bad_regions.hg19.intervals")
                                                                     
  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbsnp ROD to use (VCF)", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  val queueLogDir: String = ".qlog/"


  def script = {

    val first: Boolean = true
    val bad: Boolean = true

    val badRecalFile1: File = new File("bad_recal1.csv");
    val badRecalFile2: File = new File("bad_recal2.csv");
    val badBam: File        = new File("bad.bam");
    val badPath1: String    = "bad1";
    val badPath2: String    = "bad2";
    
    val goodRecalFile1: File = new File("good_recal1.csv")
    val goodRecalFile2: File = new File("good_recal2.csv")
    val goodBam: File        = new File("good.bam")
    val goodPath1: String    = "good1"
    val goodPath2: String    = "good2"
    
    add(cov(input, badRecalFile1, first, bad),
        recal(input, badRecalFile1, badBam, bad),
        cov(badBam, badRecalFile2, !first, bad),
        analyzeCovariates(badRecalFile1, badPath1),
        analyzeCovariates(badRecalFile2, badPath2))

    add(cov(input, goodRecalFile1, first, !bad),
        recal(input, goodRecalFile1, goodBam, !bad),
        cov(goodBam, goodRecalFile2, !first, !bad),
        analyzeCovariates(goodRecalFile1, goodPath1),
        analyzeCovariates(goodRecalFile2, goodPath2))
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = GATKjar
    this.reference_sequence = reference
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  case class cov (inBam: File, outRecalFile: File, FIRST: Boolean, BAD: Boolean) extends CountCovariates with CommandLineGATKArgs {
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    this.useOriginalQualities = FIRST
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
    if (BAD) {
      this.intervals :+= badInterval
      this.scatterCount = 85
    }
    else
      this.excludeIntervals :+= badInterval
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File, BAD: Boolean) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.out = outBam
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
    if (BAD) {
      this.intervals :+= badInterval
      this.scatterCount = 85
    }
    else
      this.excludeIntervals :+= badInterval
  }

  case class analyzeCovariates (inRecalFile: File, outPath: String) extends AnalyzeCovariates {
    this.jarFile = ACJar
    this.resources = R
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }
}