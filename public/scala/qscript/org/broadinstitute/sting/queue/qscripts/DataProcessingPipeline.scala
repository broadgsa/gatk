package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils

class DataProcessingPipeline extends QScript {
  qscript =>

  /****************************************************************************
  * Required Parameters (if default values are not good for you)
  ****************************************************************************/


  @Input(doc="input BAM file - or list of BAM files", fullName="input", shortName="i", required=true)
  var input: File = _

  @Input(doc="path to R resources folder inside the Sting repository", fullName="path_to_r", shortName="r", required=true)
  var R: String = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: File = _

  /****************************************************************************
  * Optional Parameters
  ****************************************************************************/


//  @Input(doc="path to Picard's SortSam.jar (if re-aligning a previously processed BAM file)", fullName="path_to_sort_jar", shortName="sort", required=false)
//  var sortSamJar: File = _
//

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: File = _

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Input(doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", fullName="project", shortName="p", required=false)
  var projectName: String = "project"

  @Input(doc="Output path for the processed BAM files.", fullName="output_directory", shortName="outputDir", required=false)
  var outputDir: String = ""

  @Input(doc="the -L interval string to be used by GATK - output bams at interval only", fullName="gatk_interval_string", shortName="L", required=false)
  var intervalString: String = ""

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _

  @Input(doc="Cleaning model: KNOWNS_ONLY, USE_READS or USE_SW", fullName="clean_model", shortName="cm", required=false)
  var cleaningModel: String = "USE_READS"

  @Input(doc="Decompose input BAM file and fully realign it using BWA and assume Single Ended reads", fullName="use_bwa_single_ended", shortName="bwase", required=false)
  var useBWAse: Boolean = false

  @Input(doc="Decompose input BAM file and fully realign it using BWA and assume Pair Ended reads", fullName="use_bwa_pair_ended", shortName="bwape", required=false)
  var useBWApe: Boolean = false

  @Input(doc="Number of threads BWA should use", fullName="bwa_threads", shortName="bt", required=false)
  var bwaThreads: Int = 1

  @Input(doc="Dont perform validation on the BAM files", fullName="no_validation", shortName="nv", required=false)
  var noValidation: Boolean = false


  /****************************************************************************
  * Global Variables
  ****************************************************************************/

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output
  var nContigs: Int = 0               // Use the number of contigs for scatter gathering jobs
  var cleanModelEnum: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS

  if (cleaningModel == "KNOWNS_ONLY")  {
    cleanModelEnum = ConsensusDeterminationModel.KNOWNS_ONLY
  }
  else if (cleaningModel == "USE_SW") {
    cleanModelEnum = ConsensusDeterminationModel.USE_SW
  }



  /****************************************************************************
  * Helper classes and methods
  ****************************************************************************/

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String)
  {}

  // Utility function to merge all bam files of similar samples. Generates one BAM file per sample.
  // It uses the sample information on the header of the input BAM files.
  //
  // Because the realignment only happens after these scripts are executed, in case you are using
  // bwa realignment, this function will operate over the original bam files and output over the
  // (to be realigned) bam files.
  def createSampleFiles(bamFiles: List[File], realignedBamFiles: List[File]): Map[String, File] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, List[File]]
    val realignedIterator = realignedBamFiles.iterator
    for (bam <- bamFiles) {
      val rBam = realignedIterator.next  // advance to next element in the realignedBam list so they're in sync.

      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
          sampleTable(sample) = List(rBam)
        else if ( !sampleTable(sample).contains(rBam))
          sampleTable(sample) :+= rBam
      }
    }

    println("\n\n*** INPUT FILES ***\n")
    // Creating one file for each sample in the dataset
    val sampleBamFiles = scala.collection.mutable.Map.empty[String, File]
    for ((sample, flist) <- sampleTable) {

      println(sample + ":")
      for (f <- flist)
        println (f)
      println()

      val sampleFileName = new File(qscript.outputDir + qscript.projectName + "." + sample + ".bam")
      sampleBamFiles(sample) = sampleFileName
      add(joinBams(flist, sampleFileName))
    }
    println("*** INPUT FILES ***\n\n")

    return sampleBamFiles.toMap
  }

  // Rebuilds the Read Group string to give BWA
  def addReadGroups(inBam: File, outBam: File, samReader: SAMFileReader) {
    val readGroups = samReader.getFileHeader.getReadGroups
    var index: Int = readGroups.length
    for (rg <- readGroups) {
      val intermediateInBam: File = if (index == readGroups.length) { inBam } else { swapExt(outBam, ".bam", index+1 + "-rg.bam") }
      val intermediateOutBam: File = if (index > 1) {swapExt(outBam, ".bam", index + "-rg.bam") } else { outBam}
      val readGroup = new ReadGroup(rg.getReadGroupId, rg.getPlatform, rg.getLibrary, rg.getPlatformUnit, rg.getSample, rg.getSequencingCenter, rg.getDescription)
      add(addReadGroup(intermediateInBam, intermediateOutBam, readGroup))
      index = index - 1
    }
  }

  // Takes a list of processed BAM files and realign them using the BWA option requested  (bwase or bwape).
  // Returns a list of realigned BAM files.
  def performAlignment(bams: List[File]): List[File] = {
    var realignedBams: List[File] = List()
    var index = 1
    for (bam <- bams) {
      val readSortedBam = swapExt(bam, ".bam", "." + index + ".sorted.bam" )
      val saiFile1 = swapExt(bam, ".bam", "." + index + ".1.sai")
      val saiFile2 = swapExt(bam, ".bam", "." + index + ".2.sai")
      val realignedSamFile = swapExt(bam, ".bam", "." + index + ".realigned.sam")
      val realignedBamFile = swapExt(bam, ".bam", "." + index + ".realigned.bam")
      val rgRealignedBamFile = swapExt(bam, ".bam", "." + index + ".realigned.rg.bam")
      if (useBWAse) {
        add(bwa_aln_se(bam, saiFile1),
            bwa_sam_se(bam, saiFile1, realignedSamFile))
      }
      else {
        add(sortSam(bam, readSortedBam, SortOrder.queryname),
            bwa_aln_pe(readSortedBam, saiFile1, 1),
            bwa_aln_pe(readSortedBam, saiFile2, 2),
            bwa_sam_pe(readSortedBam, saiFile1, saiFile2, realignedSamFile))
      }
      add(sortSam(realignedSamFile, realignedBamFile, SortOrder.coordinate))
      addReadGroups(realignedBamFile, rgRealignedBamFile, new SAMFileReader(bam))
      realignedBams :+= rgRealignedBamFile
      index = index + 1
    }
    return realignedBams
  }


  /****************************************************************************
  * Main script
  ****************************************************************************/


  def script = {

    // keep a record of the number of contigs in the first bam file in the list
    val bams = QScriptUtils.createListFromFile(input)
    nContigs = QScriptUtils.getNumberOfContigs(bams(0))

    val realignedBams = if (useBWApe || useBWAse) {performAlignment(bams)} else {bams}

    // Generate a BAM file per sample joining all per lane files if necessary
    val sampleBamFiles: Map[String, File] = createSampleFiles(bams, realignedBams)

    // Final output list of processed bam files
    var cohortList: List[File] = List()

    // Simple progress report
    println("\nFound the following samples: ")
    for ((sample, file) <- sampleBamFiles)
      println("\t" + sample + " -> " + file)
    println("\n")

    // If this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + projectName + ".intervals")
    if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY)
      add(target(null, globalIntervals))

    // Put each sample through the pipeline
    for ((sample, bam) <- sampleBamFiles) {

      // BAM files generated by the pipeline
      val cleanedBam = swapExt(bam, ".bam", ".clean.bam")
      val dedupedBam = swapExt(bam, ".bam", ".clean.dedup.bam")
      val recalBam   = swapExt(bam, ".bam", ".clean.dedup.recal.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val metricsFile     = swapExt(bam, ".bam", ".metrics")
      val preRecalFile    = swapExt(bam, ".bam", ".pre_recal.csv")
      val postRecalFile   = swapExt(bam, ".bam", ".post_recal.csv")
      val preOutPath      = swapExt(bam, ".bam", ".pre")
      val postOutPath     = swapExt(bam, ".bam", ".post")
      val preValidateLog  = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      if (!noValidation) {
        add(validate(bam, preValidateLog),
            validate(recalBam, postValidateLog))
      }

      if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
        add(target(bam, targetIntervals))

      add(clean(bam, targetIntervals, cleanedBam),
          dedup(cleanedBam, dedupedBam, metricsFile),
          cov(dedupedBam, preRecalFile),
          recal(dedupedBam, preRecalFile, recalBam),
          cov(recalBam, postRecalFile),
          analyzeCovariates(preRecalFile, preOutPath),
          analyzeCovariates(postRecalFile, postOutPath))


      cohortList :+= recalBam
    }

    // output a BAM list with all the processed per sample files
    val cohortFile = new File(qscript.outputDir + qscript.projectName + ".cohort.list")
    add(writeList(cohortList, cohortFile))
  }



  /****************************************************************************
  * Classes (GATK Walkers)
  ****************************************************************************/


  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }


  case class target (inBams: File, outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file :+= inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    if (!indels.isEmpty)
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
    if (!qscript.indels.isEmpty)
      this.rodBind :+= RodBind("indels", "VCF", qscript.indels)
    this.consensusDeterminationModel =  consensusDeterminationModel
    this.compress = 0
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }



  /****************************************************************************
   * Classes (non-GATK programs)
   ****************************************************************************/


  case class analyzeCovariates (inRecalFile: File, outPath: File) extends AnalyzeCovariates {
    this.resources = qscript.R
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }

  case class dedup (inBam: File, outBam: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inBam)
    this.output = outBam
    this.metrics = metricsFile
    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
  }

  case class joinBams (inBams: List[File], outBam: File) extends MergeSamFiles with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = inBams
    this.output = outBam
    this.analysisName = queueLogDir + outBam + ".joinBams"
    this.jobName = queueLogDir + outBam + ".joinBams"
  }

  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inSam)
    this.output = outBam
    this.sortOrder = sortOrderP
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input = List(inBam)
    this.output = outLog
    this.maxRecordsInRam = 100000
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = queueLogDir + outLog + ".validate"
    this.jobName = queueLogDir + outLog + ".validate"
  }


  case class addReadGroup (inBam: File, outBam: File, readGroup: ReadGroup) extends AddOrReplaceReadGroups with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inBam)
    this.output = outBam
    this.RGID = readGroup.id
    this.RGCN = readGroup.cn
    this.RGDS = readGroup.ds
    this.RGLB = readGroup.lb
    this.RGPL = readGroup.pl
    this.RGPU = readGroup.pu
    this.RGSM = readGroup.sm
    this.memoryLimit = 4
    this.isIntermediate = true
    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
  }

  case class bwa_aln_se (inBam: File, outSai: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file") var sai = outSai
    def commandLine = bwaPath + " aln -t " + bwaThreads + " -q 5 " + reference + " -b " + bam + " > " + sai
    this.analysisName = queueLogDir + outSai + ".bwa_aln_se"
    this.jobName = queueLogDir + outSai + ".bwa_aln_se"
  }

  case class bwa_aln_pe (inBam: File, outSai1: File, index: Int) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 1st mating pair") var sai = outSai1
    def commandLine = bwaPath + " aln -t " + bwaThreads + " -q 5 " + reference + " -b" + index + " " + bam + " > " + sai
    this.analysisName = queueLogDir + outSai1 + ".bwa_aln_pe1"
    this.jobName = queueLogDir + outSai1 + ".bwa_aln_pe1"
  }

  case class bwa_sam_se (inBam: File, inSai: File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file") var sai = inSai
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " samse " + reference + " " + sai + " " + bam + " > " + alignedBam
    this.analysisName = queueLogDir + outBam + ".bwa_sam_se"
    this.jobName = queueLogDir + outBam + ".bwa_sam_se"
  }

  case class bwa_sam_pe (inBam: File, inSai1: File, inSai2:File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + bam + " " + bam + " > " + alignedBam
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
