package core

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import net.sf.samtools.SAMFileReader

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class RecalibrateBaseQualities extends QScript {

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=false)
  var R: String = new File("/humgen/gsa-scr1/carneiro/stable/R")

  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbsnp ROD to use (VCF)", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  val queueLogDir: String = ".qlog/"
  var nContigs: Int = 0

  def getNumberOfContigs(bamFile: File): Int = {
    val samReader = new SAMFileReader(new File(bamFile))
    return samReader.getFileHeader.getSequenceDictionary.getSequences.size()
  }

  def script = {

    nContigs = getNumberOfContigs(input)

    val recalFile1: File = new File("recal1.csv")
    val recalFile2: File = new File("recal2.csv")
    val recalBam: File        = swapExt(input, ".bam", "recal.bam")
    val path1: String    = "before"
    val path2: String    = "after"
    
    add(cov(input, recalFile1),
        recal(input, recalFile1, recalBam),
        cov(recalBam, recalFile2),
        analyzeCovariates(recalFile1, path1),
        analyzeCovariates(recalFile2, path2))
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = GATKjar
    this.reference_sequence = reference
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
    this.scatterCount = nContigs
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.out = outBam
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
    this.scatterCount = nContigs
  }

  case class analyzeCovariates (inRecalFile: File, outPath: String) extends AnalyzeCovariates {
    this.resources = R
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }
}