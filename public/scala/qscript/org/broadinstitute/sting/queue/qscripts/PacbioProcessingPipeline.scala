package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.queue.extensions.picard.{SortSam, AddOrReplaceReadGroups}
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class RecalibrateBaseQualities extends QScript {

  @Input(doc="input FASTA/FASTQ/BAM file - or list of FASTA/FASTQ/BAM files. ", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=true)
  var R: String = _

  @Input(doc="Reference fasta file", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp VCF file to use ", shortName="D", required=true)
  var dbSNP: File = _

  @Input(doc="Number of jobs to scatter/gather. Default: 0." , shortName = "sg", required=false)
  var threads: Int = 0

  @Input(doc="Sample Name to fill in the Read Group information (only necessary if using fasta/fastq)" , shortName = "sn", required=false)
  var sample: String = "NA"

  @Input(doc="The path to the binary of bwa to align fasta/fastq files", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Hidden
  @Input(doc="The default base qualities to use before recalibration. Default is Q20 (should be good for every dataset)." , shortName = "dbq", required=false)
  var dbq: Int = 20





  val queueLogDir: String = ".qlog/"

  def script = {

    val fileList: List[File] = QScriptUtils.createListFromFile(input)

    for (file: File <- fileList) {

      var USE_BWA: Boolean = false

      println("DEBUG: processing " + file + "\nDEBUG: name -- " + file.getName)

      if (file.endsWith(".fasta") || file.endsWith(".fq")) {
        if (bwaPath == null) {
          throw new UserException("You provided a fasta/fastq file but didn't provide the path for BWA");
        }
        USE_BWA = true
      }

      // FASTA -> BAM steps
      val alignedSam: File = file.getName + ".aligned.sam"
      val sortedBam: File = swapExt(alignedSam, ".sam", ".bam")
      val qualBam: File = swapExt(sortedBam, ".bam", ".q.bam")
      val rgBam: File = swapExt(file, ".bam", ".rg.bam")

      val bamBase = if (USE_BWA) {rgBam} else {file}

      // BAM Steps
      val recalFile1: File = swapExt(bamBase, ".bam", ".recal1.csv")
      val recalFile2: File = swapExt(bamBase, ".bam", ".recal2.csv")
      val recalBam: File   = swapExt(bamBase, ".bam", ".recal.bam")
      val path1: String    = recalBam + ".before"
      val path2: String    = recalBam + ".after"


      if (USE_BWA) {
        add(align(file, alignedSam),
            sortSam(alignedSam, sortedBam),
            addQuals(sortedBam, qualBam, dbq),
            addReadGroup(qualBam, rgBam, sample))
      }

      add(cov(bamBase, recalFile1),
          recal(bamBase, recalFile1, recalBam),
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

  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = reference
  }


  case class align(@Input inFastq: File, @Output outSam: File) extends ExternalCommonArgs {
    def commandLine = bwaPath + " bwasw -b5 -q2 -r1 -z10 -t8 " + reference + " " + inFastq + " > " + outSam
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

  case class addReadGroup (inBam: File, outBam: File, sample: String) extends AddOrReplaceReadGroups with ExternalCommonArgs {
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
    this.scatterCount = threads
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.out = outBam
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
    this.scatterCount = threads
  }

  case class analyzeCovariates (inRecalFile: File, outPath: String) extends AnalyzeCovariates {
    this.resources = R
    this.recal_file = inRecalFile
    this.output_dir = outPath
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }
}
