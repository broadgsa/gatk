package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ReorderSam, SortSam, AddOrReplaceReadGroups}
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class PacbioProcessingPipeline extends QScript {

  @Input(doc="input FASTA/FASTQ/BAM file - or list of FASTA/FASTQ/BAM files. ", shortName="i", required=true)
  var input: File = _

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

  @Input(doc="Input is a BLASR generated BAM file", shortName = "blasr", fullName="blasr_bam", required=false)
  var BLASR_BAM: Boolean = false

  @Hidden
  @Input(doc="The default base qualities to use before recalibration. Default is Q20 (should be good for every dataset)." , shortName = "dbq", required=false)
  var dbq: Int = 20

  @Hidden
  @Input(shortName="bwastring", required=false)
  var bwastring: String = ""

  @Hidden
  @Input(shortName = "test", fullName = "test_mode", required = false)
  var testMode: Boolean = false

  val queueLogDir: String = ".qlog/"

  def script() {

    val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)

    for (file: File <- fileList) {

      var USE_BWA: Boolean = false
      var resetQuals: Boolean = true

      if (file.endsWith(".fasta") || file.endsWith(".fq") || file.endsWith(".fastq")) {
        if (bwaPath == null) {
          throw new UserException("You provided a fasta/fastq file but didn't provide the path for BWA");
        }
        USE_BWA = true
        if (file.endsWith(".fq") || file.endsWith(".fastq"))
          resetQuals = false
      }

      // FASTA -> BAM steps
      val alignedSam: File = file.getName + ".aligned.sam"
      val sortedBam: File = swapExt(alignedSam, ".sam", ".bam")
      val rgBam: File = swapExt(file, ".bam", ".rg.bam")

      val bamBase = if (USE_BWA) {rgBam} else {file}

      // BAM Steps
      val mqBAM: File = swapExt(bamBase, ".bam", ".mq.bam")
      val recalFile1: File = swapExt(bamBase, ".bam", ".recal1.table")
      val recalFile2: File = swapExt(bamBase, ".bam", ".recal2.table")
      val recalBam: File   = swapExt(bamBase, ".bam", ".recal.bam")
      val path1: String    = recalBam + ".before"
      val path2: String    = recalBam + ".after"

      if (USE_BWA) {
        add(align(file, alignedSam),
            sortSam(alignedSam, sortedBam),
            addReadGroup(sortedBam, rgBam, sample))
      }

      else if (BLASR_BAM) {
        val reorderedBAM = swapExt(bamBase, ".bam", ".reordered.bam")
        add(reorder(bamBase, reorderedBAM),
            changeMQ(reorderedBAM, mqBAM))
      }

      val bam = if (BLASR_BAM) {mqBAM} else {bamBase}

      add(cov(bam, recalFile1, resetQuals),
          recal(bam, recalFile1, recalBam),
          cov(recalBam, recalFile2, false))
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
    def commandLine = bwaPath + " bwasw -b5 -q2 -r1 -z20 -t16 " + reference + " " + inFastq + " > " + outSam
    this.memoryLimit = 8
    this.analysisName = queueLogDir + outSam + ".bwa_sam_se"
    this.jobName = queueLogDir + outSam + ".bwa_sam_se"
  }

  case class sortSam (inSam: File, outBam: File) extends SortSam with ExternalCommonArgs {
    this.input = List(inSam)
    this.output = outBam
    this.memoryLimit = 8
    this.sortOrder = SortOrder.coordinate
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class reorder (inSam: File, outSam: File) extends ReorderSam with ExternalCommonArgs {
    this.input = List(inSam)
    this.output = outSam
    this.sortReference = reference
  }

  case class changeMQ(inBam: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outBam
    this.read_filter :+= "ReassignMappingQuality"
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

  case class cov (inBam: File, outRecalFile: File, resetQuals: Boolean) extends BaseQualityScoreRecalibrator with CommandLineGATKArgs {
    if (resetQuals) 
      this.DBQ = dbq
    this.knownSites :+= dbSNP
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.out = outRecalFile
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
    this.scatterCount = threads
    this.read_filter :+= "BadCigar"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.DBQ = dbq
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.out = outBam
    this.disable_indel_quals = true
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
    this.read_filter :+= "BadCigar"
    this.scatterCount = threads
  }
}
