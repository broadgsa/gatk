package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.queue.extensions.picard.{SortSam, AddOrReplaceReadGroups}

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class RecalibrateBaseQualities extends QScript {

  @Input(doc="input FASTA file, BAM file - or list of FASTA/BAM files. ", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=true)
  var R: String = _

  @Input(doc="Reference fasta file", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp VCF file to use ", shortName="D", required=true)
  var dbSNP: File = _

  @Input(doc="Default base qualities. Overrides the file's original base qualities with given value. Must be used if the file does not have base qualities." , shortName = "dbq", required=false)
  var dbq: Int = -1

  @Input(doc="Number of jobs to scatter/gather. Default is the number of contigs in the dataset" , shortName = "sg", required=false)
  var threads: Int = -1

  @Input(doc="Sample Name" , shortName = "sn", required=false)
  var sample: String = ""

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _





  val queueLogDir: String = ".qlog/"
  var nContigs: Int = 0
  var ADD_BASE_QUALITIES = false

  def script = {

    if (dbq >= 0)
      ADD_BASE_QUALITIES = true

    val fileList = QScriptUtils.createListFromFile(input)
    nContigs = if (threads >= 0) {threads} else {QScriptUtils.getNumberOfContigs(fileList(0))}

    for (file <- fileList) {
      val qualBam: File = swapExt(file, ".bam", ".quals.bam")
      val rgBam: File = if (ADD_BASE_QUALITIES)  {swapExt(file, ".bam", ".rg.bam")} else {file}
      val recalFile1: File = swapExt(file, ".bam", ".recal1.csv")
      val recalFile2: File = swapExt(file, ".bam", ".recal2.csv")
      val recalBam: File   = swapExt(file, ".bam", ".recal.bam")
      val path1: String    = recalBam + ".before"
      val path2: String    = recalBam + ".after"


      if (ADD_BASE_QUALITIES) {
        add(addQuals(file, qualBam, dbq),
            addReadGroup(qualBam, rgBam, sample))
      }

      add(cov(rgBam, recalFile1),
          recal(rgBam, recalFile1, recalBam),
          cov(recalBam, recalFile2),
          analyzeCovariates(recalFile1, path1),
          analyzeCovariates(recalFile2, path2))
    }
  }


  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
  }


  case class align(@Input inFastq: File, @Output outSam: File) extends ExternalCommonArgs {
    def commandLine = bwaPath + " bwasw " + reference + " " + inFastq + " > " + outSam
    this.analysisName = queueLogDir + outSam + ".bwa_sam_se"
    this.jobName = queueLogDir + outSam + ".bwa_sam_se"
  }

  case class sortSam (@Input inSam: File, @Output outBam: File) extends SortSam with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inSam)
    this.output = outBam
    this.sortOrder = SortOrder.coordinate
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class addQuals(inBam: File, outBam: File, qual: Int) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outBam
    this.DBQ = qual
  }

  case class addReadGroup (inBam: File, outBam: File, sample: String) extends AddOrReplaceReadGroups {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inBam)
    this.output = outBam
    this.RGID = "1"
    this.RGCN = "BI"
    this.RGPL = "PacBio_RS"
    this.RGSM = sample
    this.RGLB = "default_library"
    this.RGPU = "default_pu"
    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
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
    this.output_dir = outPath
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }
}
