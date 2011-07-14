package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import net.sf.samtools.SAMFileReader
import io.Source._
import org.broadinstitute.sting.queue.qscripts.utils.Utils

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

  def script = {

    val bamList = Utils.createListFromFile(input)
    nContigs = Utils.getNumberOfContigs(bamList(0))

    for (bam <- bamList) {

      val recalFile1: File = swapExt(bam, ".bam", ".recal1.csv")
      val recalFile2: File = swapExt(bam, ".bam", ".recal2.csv")
      val recalBam: File   = swapExt(bam, ".bam", ".recal.bam")
      val path1: String    = bam + "before"
      val path2: String    = bam + "after"

      add(cov(bam, recalFile1),
          recal(bam, recalFile1, recalBam),
          cov(recalBam, recalFile2),
          analyzeCovariates(recalFile1, path1),
          analyzeCovariates(recalFile2, path2))
    }
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
