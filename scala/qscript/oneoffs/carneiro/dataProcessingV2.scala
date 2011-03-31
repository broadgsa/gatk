package oneoffs.carneiro

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.ListWriterFunction

import net.sf.samtools.{SAMFileReader,SAMReadGroupRecord}

import scala.io.Source._
import collection.JavaConversions._


class dataProcessingV2 extends QScript {
  qscript =>

  /****************************************************************************
  * Required Parameters (if default values are not good for you)
  ****************************************************************************/


  @Input(doc="input BAM file - or list of BAM files", fullName="input", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to GenomeAnalysisTK.jar", fullName="path_to_gatk_jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="path to AnalyzeCovariates.jar", fullName="path_to_ac_jar", shortName="ac", required=true)
  var ACJar: File = _

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=true)
  var R: String = _

  @Input(doc="path to Picard's MarkDuplicates.jar", fullName="path_to_dedup_jar", shortName="dedup", required=false)
  var dedupJar: File = new File("/seq/software/picard/current/bin/MarkDuplicates.jar")

  @Input(doc="path to Picard's MergeSamFiles.jar", fullName="path_to_merge_jar", shortName="merge", required=false)
  var mergeBamJar: File = new File("/seq/software/picard/current/bin/MergeSamFiles.jar")

  @Input(doc="path to Picard's ValidateSamFile.jar", fullName="path_to_validate_jar", shortName="validate", required=false)
  var validateSamJar: File = new File("/seq/software/picard/current/bin/ValidateSamFile.jar")

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")



  /****************************************************************************
  * Optional Parameters
  ****************************************************************************/


  @Input(doc="path to Picard's RevertSam.jar (if re-aligning a previously processed BAM file)", fullName="path_to_revert_jar", shortName="revert", required=false)
  var revertSamJar: File = _

  @Input(doc="path to Picard's SortSam.jar (if re-aligning a previously processed BAM file)", fullName="path_to_sort_jar", shortName="sort", required=false)
  var sortSamJar: File = _

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf")

  @Input(doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", fullName="project", shortName="p", required=false)
  var projectName: String = "project"

  @Input(doc="Output path for the processed BAM files.", fullName="output_directory", shortName="outputDir", required=false)
  var outputDir: String = ""

  @Input(doc="the -L interval string to be used by GATK - output bams at interval only", fullName="gatk_interval_string", shortName="L", required=false)
  var intervalString: String = ""

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _

  @Input(doc="Perform cleaning on knowns only", fullName="knowns_only", shortName="knowns", required=false)
  var knownsOnly: Boolean = false

  @Input(doc="Perform cleaning using Smith Waterman", fullName="use_smith_waterman", shortName="sw", required=false)
  var useSW: Boolean = false

  @Input(doc="Decompose input BAM file and fully realign it using BWA and assume Single Ended reads", fullName="use_bwa_single_ended", shortName="bwase", required=false)
  var useBWAse: Boolean = false

  @Input(doc="Decompose input BAM file and fully realign it using BWA and assume Pair Ended reads", fullName="use_bwa_pair_ended", shortName="bwape", required=false)
  var useBWApe: Boolean = false



  /****************************************************************************
  * Global Variables
  ****************************************************************************/

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output
  var nContigs: Int = 0               // Use the number of contigs for scatter gathering jobs



  /****************************************************************************
  * Helper functions
  ****************************************************************************/


  // Utility function to check if there are multiple samples in a BAM file (currently we can't deal with that)
  def hasMultipleSamples(readGroups: java.util.List[SAMReadGroupRecord]): Boolean = {
    var sample: String = ""
    for (r <- readGroups) {
      if (sample.isEmpty)
        sample = r.getSample()
      else if (sample != r.getSample())
          return true;
    }
    return false
  }

  // Utility function to merge all bam files of similar samples. Generates on BAM file per sample.
  // It uses the sample information on the header of the input BAM files.
  def createSampleFiles(bamFiles: List[File]): Map[String, File] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, List[File]]
    for (bam <- bamFiles) {
      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader()
      val readGroups = header.getReadGroups()

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample()
        if (!sampleTable.contains(sample))
          sampleTable(sample) = List(bam)
        else if ( !sampleTable(sample).contains(bam))
          sampleTable(sample) :+= bam
      }
    }

    // Creating one file for each sample in the dataset
    val sampleBamFiles = scala.collection.mutable.Map.empty[String, File]
    for ((sample, flist) <- sampleTable) {
      val sampleFileName = new File(qscript.outputDir + qscript.projectName + "." + sample + ".bam")
      sampleBamFiles(sample) = sampleFileName
      add(joinBams(flist, sampleFileName))
    }
    return sampleBamFiles.toMap
  }

  // Checks how many contigs are in the dataset. Uses the BAM file header information.
  def getNumberOfContigs(): Int = {
    val bam = fromFile(input).getLines.next
    val samReader = new SAMFileReader(new File(bam))
    return samReader.getFileHeader.getSequenceDictionary.getSequences.size()
  }


  // Takes a list of processed BAM files, revert them to unprocessed and realigns each lane, producing a list of
  // per-lane aligned bam files, ready to be processed.
  def performAlignment(bamList: File): List[File] = {
    return List()
  }

  def createListFromFile(in: File):List[File] = {
    if (in.toString.endsWith("bam"))
      return List(in)
    var l: List[File] = List()
    for (bam <- fromFile(in).getLines)
      l :+= new File(bam)
    return l
  }



  /****************************************************************************
  * Main script
  ****************************************************************************/


  def script = {

    //todo -- (option - BWA) run BWA on each bam file (per lane bam file) before performing per sample processing
    val perLaneAlignedBamFiles: List[File] = if (useBWApe || useBWAse) { performAlignment(input) } else { createListFromFile(input) }

    // Generate a BAM file per sample joining all per lane files if necessary
    val sampleBamFiles = createSampleFiles(perLaneAlignedBamFiles)

    // keep a record of the number of contigs in this bam file (they should all be the same
    nContigs = getNumberOfContigs()

    // Final output list of processed bam files
    var cohortList: List[File] = List()

    // Simple progress report
    println("\nFound the following samples: ")
    for ((sample, file) <- sampleBamFiles)
      println("\t" + sample + " -> " + file)

    // If this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + projectName + ".intervals")
    if (knownsOnly)
      add(target(null, globalIntervals))

    // Put each sample through the pipeline
    for ((sample, bam) <- sampleBamFiles) {

      // BAM files generated by the pipeline
      val cleanedBam = swapExt(bam, ".bam", ".clean.bam")
      val dedupedBam = swapExt(bam, ".bam", ".clean.dedup.bam")
      val recalBam   = swapExt(bam, ".bam", ".clean.dedup.recal.bam")

      // Accessory files
      val targetIntervals = if (knownsOnly) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val metricsFile     = swapExt(bam, ".bam", ".metrics")
      val preRecalFile    = swapExt(bam, ".bam", ".pre_recal.csv")
      val postRecalFile   = swapExt(bam, ".bam", ".post_recal.csv")
      val preOutPath      = swapExt(bam, ".bam", ".pre")
      val postOutPath     = swapExt(bam, ".bam", ".post")
      val preValidateLog  = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      add(validate(bam, preValidateLog))

      if (!knownsOnly)
        add(target(bam, targetIntervals))

      add(clean(bam, targetIntervals, cleanedBam),
          dedup(cleanedBam, dedupedBam, metricsFile),
          cov(dedupedBam, preRecalFile),
          recal(dedupedBam, preRecalFile, recalBam),
          cov(recalBam, postRecalFile),
          analyzeCovariates(preRecalFile, preOutPath),
          analyzeCovariates(postRecalFile, postOutPath),
          validate(recalBam, postValidateLog))

      cohortList :+= recalBam
    }

    // output a BAM list with all the processed per sample files
    val cohortFile = new File(qscript.outputDir + qscript.projectName + ".cohort.list")
    add(writeList(cohortList, cohortFile))
  }



  /****************************************************************************
  * Classes (Walkers and non-GATK programs)
  ****************************************************************************/


  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.GATKjar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  case class target (inBams: File, outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (!knownsOnly)
      this.input_file :+= inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.rodBind :+= RodBind("indels", "VCF", indels)
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class clean (inBams: File, tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file :+= inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.rodBind :+= RodBind("indels", "VCF", qscript.indels)
    this.useOnlyKnownIndels = knownsOnly
    this.doNotUseSW = useSW
    this.compress = 0
    this.U = org.broadinstitute.sting.gatk.arguments.ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION  // todo -- update this with the last consensus between Tim, Matt and Eric. This is ugly!
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  //todo -- add scatter gather capability (waiting for khalid's modifications to the queue base
  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    this.scatterCount = nContigs
    this.jobQueue = "gsa"      // should take this out once scatter gather is available.
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  //todo -- add scatter gather capability (waiting for khalid's modifications to the queue base
  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    @Output(doc="recalibrated bam index") var recalIndex = new File(outBam + ".bai")
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.U = org.broadinstitute.sting.gatk.arguments.ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION  // todo -- update this with the last consensus between Tim, Matt and Eric. This is ugly!
    this.index_output_bam_on_the_fly = true // todo -- implemente @gather for BAM index
//    this.scatterCount = nContigs
    this.isIntermediate = false
    this.jobQueue = "gsa"      // should take this out once scatter gather is available.
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"

  }



  // Outside tools (not GATK walkers)

  case class analyzeCovariates (inRecalFile: File, outPath: File) extends AnalyzeCovariates {
    this.jarFile = qscript.ACJar
    this.resources = qscript.R
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }

  case class dedup (inBam: File, outBam: File, metricsFile: File) extends PicardBamFunction {
    @Input(doc="fixed bam") var clean = inBam
    @Output(doc="deduped bam") var deduped = outBam
    @Output(doc="deduped bam index") var dedupedIndex = new File(outBam + "bai")
    @Output(doc="metrics file") var metrics = metricsFile
    override def inputBams = List(clean)
    override def outputBam = deduped
    override def commandLine = super.commandLine + " M=" + metricsFile
    this.sortOrder = null
    this.createIndex = true
    this.memoryLimit = 6
    this.isIntermediate = true
    this.jarFile = qscript.dedupJar
    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
  }

  case class joinBams (inBams: List[File], outBam: File) extends PicardBamFunction {
    @Input(doc="input bam list") var join = inBams
    @Output(doc="joined bam") var joined = outBam
    @Output(doc="joined bam index") var joinedIndex = new File(outBam + "bai")
    override def inputBams = join
    override def outputBam = joined
    override def commandLine = super.commandLine + " CREATE_INDEX=true"
    this.jarFile = qscript.mergeBamJar
    this.isIntermediate = true
    this.analysisName = queueLogDir + outBam + ".joinBams"
    this.jobName = queueLogDir + outBam + ".joinBams"
  }

  case class sortSam (inSam: File, outBam: File) extends PicardBamFunction {
    @Input(doc="input unsorted sam file") var sam = inSam
    @Output(doc="sorted bam") var bam = outBam
    @Output(doc="sorted bam index") var bamIndex = new File(outBam + "bai")
    override def inputBams = List(sam)
    override def outputBam = bam
    override def commandLine = super.commandLine + " CREATE_INDEX=true"
    this.jarFile = qscript.sortSamJar
    this.isIntermediate = true
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends PicardBamFunction {
    @Input(doc="input bam list") var toValidate = inBam
    @Output(doc="validation log") var validate = outLog
    override def inputBams = List(inBam)
    override def outputBam = outLog
    override def commandLine = super.commandLine + " VALIDATE_INDEX=true MODE=SUMMARY REFERENCE_SEQUENCE=" + qscript.reference
    sortOrder = null
    this.jarFile = qscript.validateSamJar
    this.isIntermediate = false
    this.analysisName = queueLogDir + outLog + ".validate"
    this.jobName = queueLogDir + outLog + ".validate"
  }

  case class revert (inBam: File, outBam: File) extends PicardBamFunction {
    @Input(doc="old annotated bam") var oldBam = inBam
    @Output(doc="reverted bam") var revertedBam = outBam
    @Output(doc="reverted bam index") var revertedBamIndex = new File(outBam + ".bai")
    override def inputBams = List(oldBam)
    override def outputBam = revertedBam
    override def commandLine = super.commandLine + " CREATE_INDEX=true"
    this.isIntermediate = true
    this.jarFile = qscript.dedupJar
    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
  }

  case class bwa_aln_se (inBam: File, outSai: File) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file") var sai = outSai
    def commandLine = bwaPath + " aln -q 5 " + reference + " -b " + bam + " > " + sai
    this.isIntermediate = true
    this.analysisName = queueLogDir + outSai + ".bwa_aln_se"
    this.jobName = queueLogDir + outSai + ".bwa_aln_se"
  }

  case class bwa_aln_pe1 (inBam: File, outSai1: File) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 1st mating pair") var sai = outSai1
    def commandLine = bwaPath + " aln -q 5 " + reference + " -b1 " + bam + " > " + sai
    this.isIntermediate = true
    this.analysisName = queueLogDir + outSai1 + ".bwa_aln_pe1"
    this.jobName = queueLogDir + outSai1 + ".bwa_aln_pe1"
  }

  case class bwa_aln_pe2 (inBam: File, outSai2: File) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 2nd mating pair") var sai = outSai2
    def commandLine = bwaPath + " aln -q 5 " + reference + " -b2 " + bam + " > " + sai
    this.isIntermediate = true
    this.analysisName = queueLogDir + outSai2 + ".bwa_aln_pe2"
    this.jobName = queueLogDir + outSai2 + ".bwa_aln_pe2"
  }

  case class bwa_sam_se (inBam: File, inSai: File, outBam: File, readGroup: String) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file") var sai = inSai
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " samse " + reference + " " + sai + " " + bam + " -r " + readGroup + " > " + alignedBam
    this.isIntermediate = true
    this.analysisName = queueLogDir + outBam + ".bwa_sam_se"
    this.jobName = queueLogDir + outBam + ".bwa_sam_se"
  }

  case class bwa_sam_pe (inBam: File, inSai1: File, inSai2:File, outBam: File, readGroup: String) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " samse " + reference + " " + sai1 + " " + sai2 + " " + bam + " -r " + readGroup + " > " + alignedBam
    this.isIntermediate = true
    this.analysisName = queueLogDir + outBam + ".bwa_sam_pe"
    this.jobName = queueLogDir + outBam + ".bwa_sam_pe"
  }

  case class writeList(inBams: List[File], outBamList: File) extends ListWriterFunction {
    this.inputFiles = inBams
    this.listFile = outBamList
    this.analysisName = queueLogDir + outBamList + ".bamList"
    this.jobName = queueLogDir + outBamList + ".bamList"
  }


}
