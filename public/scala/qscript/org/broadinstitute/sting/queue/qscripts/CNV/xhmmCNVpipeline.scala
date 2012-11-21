package org.broadinstitute.sting.queue.qscripts.CNV

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import org.broadinstitute.sting.queue.util.DoC._
import org.broadinstitute.sting.commandline.Hidden
import java.io.{PrintStream, PrintWriter}
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.sting.gatk.walkers.coverage.CoverageUtils

class xhmmCNVpipeline extends QScript {
  qscript =>

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Input(doc = "xhmm executable file", shortName = "xhmmExec", required = true)
  var xhmmExec: File = _

  @Input(doc = "Plink/Seq executable file", shortName = "pseqExec", required = true)
  var pseqExec: File = _

  @Argument(doc = "Plink/Seq SEQDB file (Reference genome sequence)", shortName = "SEQDB", required = true)
  var pseqSeqDB: String = _

  @Input(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Input(shortName = "L", doc = "Intervals", required = false)
  var intervals: File = _

  @Argument(doc = "level of parallelism for BAM DoC.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCountInput = 0

  @Argument(doc = "Samples to run together for DoC.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Base name for files to output", shortName = "o", required = true)
  var outputBase: File = _

  @Hidden
  @Argument(doc = "How should overlapping reads from the same fragment be handled?", shortName = "countType", required = false)
  var countType = CoverageUtils.CountPileupType.COUNT_FRAGMENTS

  @Argument(doc = "Maximum depth (before GATK down-sampling kicks in...)", shortName = "MAX_DEPTH", required = false)
  var MAX_DEPTH = 20000

  @Hidden
  @Argument(doc = "Number of read-depth bins", shortName = "NUM_BINS", required = false)
  var NUM_BINS = 200

  @Hidden
  @Argument(doc = "Starting value of read-depth bins", shortName = "START_BIN", required = false)
  var START_BIN = 1

  @Argument(doc = "Minimum read mapping quality", shortName = "MMQ", required = false)
  var minMappingQuality = 0

  @Argument(doc = "Minimum base quality to be counted in depth", shortName = "MBQ", required = false)
  var minBaseQuality = 0

  @Argument(doc = "Memory (in GB) required for storing the whole matrix in memory", shortName = "wholeMatrixMemory", required = false)
  var wholeMatrixMemory = -1

  @Argument(shortName = "minTargGC", doc = "Exclude all targets with GC content less than this value", required = false)
  var minTargGC : Double = 0.1

  @Argument(shortName = "maxTargGC", doc = "Exclude all targets with GC content greater than this value", required = false)
  var maxTargGC : Double = 0.9

  @Argument(shortName = "minTargRepeats", doc = "Exclude all targets with % of repeat-masked bases less than this value", required = false)
  var minTargRepeats : Double = 0.0

  @Argument(shortName = "maxTargRepeats", doc = "Exclude all targets with % of repeat-masked bases greater than this value", required = false)
  var maxTargRepeats : Double = 0.1

  @Argument(shortName = "sampleIDsMap", doc = "File mapping BAM sample IDs to desired sample IDs", required = false)
  var sampleIDsMap: String = ""

  @Argument(shortName = "sampleIDsMapFromColumn", doc = "Column number of OLD sample IDs to map", required = false)
  var sampleIDsMapFromColumn = 1

  @Argument(shortName = "sampleIDsMapToColumn", doc = "Column number of NEW sample IDs to map", required = false)
  var sampleIDsMapToColumn = 2

  @Argument(shortName = "rawFilters", doc = "xhmm command-line parameters to filter targets and samples from raw data", required = false)
  var targetSampleFiltersString: String = ""

  @Argument(shortName = "PCAnormalize", doc = "xhmm command-line parameters to Normalize data using PCA information", required = false)
  var PCAnormalizeMethodString: String = ""

  @Argument(shortName = "normalizedFilters", doc = "xhmm command-line parameters to filter targets and samples from PCA-normalized data", required = false)
  var targetSampleNormalizedFiltersString: String = ""

  @Argument(shortName = "xhmmParams", doc = "xhmm model parameters file", required = true)
  var xhmmParamsArg: File = _

  @Argument(shortName = "discoverParams", doc = "xhmm command-line parameters for discovery step", required = false)
  var discoverCommandLineParams: String = ""

  @Argument(shortName = "genotypeParams", doc = "xhmm command-line parameters for genotyping step", required = false)
  var genotypeCommandLineParams: String = ""

  @Argument(shortName = "genotypeSubsegments", doc = "Should we also genotype all subsegments of the discovered CNV?", required = false)
  var genotypeSubsegments: Boolean = false

  @Argument(shortName = "maxTargetsInSubsegment", doc = "If genotypeSubsegments, then only consider sub-segments consisting of this number of targets or fewer", required = false)
  var maxTargetsInSubsegment = 30

  @Argument(shortName = "subsegmentGenotypeThreshold", doc = "If genotypeSubsegments, this is the default genotype quality threshold for the sub-segments", required = false)
  var subsegmentGenotypeThreshold = 20.0

  @Argument(shortName = "longJobQueue", doc = "Job queue to run the 'long-running' commands", required = false)
  var longJobQueue: String = ""


  val PREPARED_TARGS_SUFFIX: String = ".merged.interval_list"

  val RD_OUTPUT_SUFFIX: String = ".RD.txt"

  val TARGS_GC_SUFFIX = ".locus_GC.txt"
  val EXTREME_GC_TARGS_SUFFIX = ".extreme_gc_targets.txt"

  val TARGS_REPEAT_COMPLEXITY_SUFFIX = ".locus_complexity.txt"
  val EXTREME_REPEAT_COMPLEXITY_SUFFIX = ".extreme_complexity_targets.txt"

  val FILTERED_TARGS_SUFFIX: String = ".filtered_targets.txt"
  val FILTERED_SAMPS_SUFFIX: String =  ".filtered_samples.txt"


  trait WholeMatrixMemoryLimit extends CommandLineFunction {
    // Since loading ALL of the data can take significant memory:
    if (wholeMatrixMemory < 0) {
      this.memoryLimit = 24
    }
    else {
      this.memoryLimit = wholeMatrixMemory
    }
  }

  trait LongRunTime extends CommandLineFunction {
    if (longJobQueue != "")
      this.jobQueue = longJobQueue
  }

  def script = {
    val prepTargets = new PrepareTargets(List(qscript.intervals), outputBase.getPath + PREPARED_TARGS_SUFFIX, xhmmExec, referenceFile)
    add(prepTargets)

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.intervals :+= prepTargets.out
      this.jarFile = qscript.gatkJarFile
      this.reference_sequence = qscript.referenceFile
      this.logging_level = "INFO"
    }

    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = VCF_BAM_utilities.getMapOfBAMsForSample(VCF_BAM_utilities.parseBAMsInput(bams))
    val samples: List[String] = sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val groups: List[Group] = buildDoCgroups(samples, sampleToBams, samplesPerJob, outputBase)
    var docs: List[DoC] = List[DoC]()
    for (group <- groups) {
      Console.out.printf("Group is %s%n", group)
      docs ::= new DoC(group.bams, group.DoC_output, countType, MAX_DEPTH, minMappingQuality, minBaseQuality, scatterCountInput, START_BIN, NUM_BINS, Nil) with CommandLineGATKArgs
    }
    addAll(docs)

    val mergeDepths = new MergeGATKdepths(docs.map(u => u.intervalSampleOut), outputBase.getPath + RD_OUTPUT_SUFFIX, "_mean_cvg", xhmmExec, sampleIDsMap, sampleIDsMapFromColumn, sampleIDsMapToColumn, None, false) with WholeMatrixMemoryLimit
    add(mergeDepths)

    var excludeTargets : List[File] = List[File]()
    if (minTargGC > 0 || maxTargGC < 1) {
      val calcGCcontents = new GCContentByInterval with CommandLineGATKArgs
      calcGCcontents.out = outputBase.getPath + TARGS_GC_SUFFIX
      add(calcGCcontents)

      val excludeTargetsBasedOnGC = new ExcludeTargetsBasedOnValue(calcGCcontents.out, EXTREME_GC_TARGS_SUFFIX, minTargGC, maxTargGC)
      add(excludeTargetsBasedOnGC)
      excludeTargets ::= excludeTargetsBasedOnGC.out
    }

    class CalculateRepeatComplexity(outFile : String) extends CommandLineFunction {
      @Input(doc="")
      var intervals: File = prepTargets.out

      @Output(doc="")
      var out : File = new File(outFile)

      val regFile : String = outputBase.getPath + ".targets.reg"
      val locDB : String = outputBase.getPath + ".targets.LOCDB"

      val removeFiles = "rm -f " + regFile + " " + locDB
      val createRegFile = "cat " + intervals + " | awk 'BEGIN{OFS=\"\\t\"; print \"#CHR\\tBP1\\tBP2\\tID\"} {split($1,a,\":\"); chr=a[1]; if (match(chr,\"chr\")==0) {chr=\"chr\"chr} split(a[2],b,\"-\"); bp1=b[1]; bp2=bp1; if (length(b) > 1) {bp2=b[2]} print chr,bp1,bp2,NR}' > " + regFile
      val createLOCDB = pseqExec + " . loc-load --locdb " + locDB + " --file " + regFile + " --group targets --out " + locDB + ".loc-load"
      val calcRepeatMaskedPercent = pseqExec + " . loc-stats --locdb " + locDB + " --group targets --seqdb " + pseqSeqDB + " --out " + locDB + ".loc-stats"
      val extractRepeatMaskedPercent = "cat " + locDB + ".loc-stats.locstats | awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste " + intervals + " - | awk '{print $1\"\\t\"$2}' > " + out

      var command: String =
        removeFiles +
          " && " + createRegFile +
          " && " + createLOCDB +
          " && " + calcRepeatMaskedPercent +
          " && " + extractRepeatMaskedPercent

      def commandLine = command

      override def description = "Calculate the percentage of each target that is repeat-masked in the reference sequence: " + command
    }

    if (minTargRepeats > 0 || maxTargRepeats < 1) {
      val calcRepeatComplexity = new CalculateRepeatComplexity(outputBase.getPath + TARGS_REPEAT_COMPLEXITY_SUFFIX)
      add(calcRepeatComplexity)

      val excludeTargetsBasedOnRepeats = new ExcludeTargetsBasedOnValue(calcRepeatComplexity.out, EXTREME_REPEAT_COMPLEXITY_SUFFIX, minTargRepeats, maxTargRepeats)
      add(excludeTargetsBasedOnRepeats)
      excludeTargets ::= excludeTargetsBasedOnRepeats.out
    }

    val filterCenterDepths = new FilterCenterRawMatrix(mergeDepths.mergedDoC, excludeTargets)
    add(filterCenterDepths)

    val pca = new PCA(filterCenterDepths.filteredCentered)
    add(pca)

    val normalize = new Normalize(pca)
    add(normalize)

    val filterZscore = new FilterAndZscoreNormalized(normalize.normalized)
    add(filterZscore)

    val filterOriginal = new FilterOriginalData(mergeDepths.mergedDoC, filterCenterDepths, filterZscore)
    add(filterOriginal)

    val discover = new DiscoverCNVs(filterZscore.filteredZscored, filterOriginal.sameFiltered)
    add(discover)

    val genotype = new GenotypeCNVs(filterZscore.filteredZscored, discover.xcnv, filterOriginal.sameFiltered)
    add(genotype)

    if (genotypeSubsegments) {
      val genotypeSegs = new GenotypeCNVandSubsegments(filterZscore.filteredZscored, discover.xcnv, filterOriginal.sameFiltered)
      add(genotypeSegs)
    }
  }

  class ExcludeTargetsBasedOnValue(locus_valueIn : File, outSuffix : String, minVal : Double, maxVal : Double) extends InProcessFunction {
    @Input(doc="")
    var locus_value : File = locus_valueIn

    @Output(doc="")
    var out : File = new File(outputBase.getPath + outSuffix)

    def run = {
      var outWriter = new PrintWriter(new PrintStream(out))
      var elems = asScalaIterator(new XReadLines(locus_value))

      while (elems.hasNext) {
        val line = elems.next
        val splitLine = line.split("\\s+")
        val locus = splitLine(0)
        val locValStr = splitLine(1)
        try {
          val locVal = locValStr.toDouble
          if (locVal < minVal || locVal > maxVal)
            outWriter.printf("%s%n", locus)
        }
        catch {
          case nfe: NumberFormatException => println("Ignoring non-numeric value " + locValStr + " for locus " + locus)
          case e: Exception => throw e
        }
      }

      outWriter.close
    }
  }

  class FilterCenterRawMatrix(inputParam: File, excludeTargetsIn : List[File]) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val excludeTargets = excludeTargetsIn

    @Output
    val filteredCentered: File = new File(outputBase.getPath + ".filtered_centered" + RD_OUTPUT_SUFFIX)
    @Output
    val filteredTargets: File = new File(filteredCentered.getPath + FILTERED_TARGS_SUFFIX)
    @Output
    val filteredSamples: File = new File(filteredCentered.getPath + FILTERED_SAMPS_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      " --centerData --centerType target" +
      " -o " + filteredCentered +
      " --outputExcludedTargets " + filteredTargets +
      " --outputExcludedSamples " + filteredSamples
    command += excludeTargets.map(u => " --excludeTargets " + u).reduceLeft(_ + "" + _)
    if (targetSampleFiltersString != "")
      command += " " + targetSampleFiltersString

    def commandLine = command

    override def description = "Filters samples and targets and then mean-centers the targets: " + command
  }

  class PCA(inputParam: File) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    val PCAbase: String = outputBase.getPath + ".RD_PCA"

    @Output
    val outPC: File = new File(PCAbase + ".PC.txt")
    @Output
    val outPC_SD: File = new File(PCAbase + ".PC_SD.txt")
    @Output
    val outPC_LOADINGS: File = new File(PCAbase + ".PC_LOADINGS.txt")

    var command: String =
      xhmmExec + " --PCA" +
      " -r " + input +
      " --PCAfiles " + PCAbase

    def commandLine = command

    override def description = "Runs PCA on mean-centered data: " + command
  }

  class Normalize(pca: PCA) extends CommandLineFunction {
    @Input(doc = "")
    val input = pca.input

    @Input(doc = "")
    val inPC = pca.outPC

    @Input(doc = "")
    val inPC_SD = pca.outPC_SD

    @Input(doc = "")
    val inPC_LOADINGS = pca.outPC_LOADINGS

    @Output
    val normalized: File = new File(outputBase.getPath + ".PCA_normalized.txt")

    var command: String =
      xhmmExec + " --normalize" +
      " -r " + input +
      " --PCAfiles " + pca.PCAbase +
      " --normalizeOutput " + normalized
    if (PCAnormalizeMethodString != "")
      command += " " + PCAnormalizeMethodString

    def commandLine = command

    override def description = "Normalizes mean-centered data using PCA information: " + command
  }

  class FilterAndZscoreNormalized(inputParam: File) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Output
    val filteredZscored: File = new File(outputBase.getPath + ".PCA_normalized.filtered.sample_zscores" + RD_OUTPUT_SUFFIX)
    @Output
    val filteredTargets: File = new File(filteredZscored.getPath + FILTERED_TARGS_SUFFIX)
    @Output
    val filteredSamples: File = new File(filteredZscored.getPath + FILTERED_SAMPS_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      " --centerData --centerType sample --zScoreData" +
      " -o " + filteredZscored +
      " --outputExcludedTargets " + filteredTargets +
      " --outputExcludedSamples " + filteredSamples
    if (targetSampleNormalizedFiltersString != "")
      command += " " + targetSampleNormalizedFiltersString

    def commandLine = command

    override def description = "Filters and z-score centers (by sample) the PCA-normalized data: " + command
  }

  class FilterOriginalData(inputParam: File, filt1: FilterCenterRawMatrix, filt2: FilterAndZscoreNormalized) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val targFilters: List[File] = List(filt1.filteredTargets, filt2.filteredTargets).map(u => new File(u))

    @Input(doc = "")
    val sampFilters: List[File] = List(filt1.filteredSamples, filt2.filteredSamples).map(u => new File(u))

    @Output
    val sameFiltered: File = new File(outputBase.getPath + ".same_filtered" + RD_OUTPUT_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      targFilters.map(u => " --excludeTargets " + u).reduceLeft(_ + "" + _) +
      sampFilters.map(u => " --excludeSamples " + u).reduceLeft(_ + "" + _) +
      " -o " + sameFiltered

    def commandLine = command

    override def description = "Filters original read-depth data to be the same as filtered, normalized data: " + command
  }

  class DiscoverCNVs(inputParam: File, origRDParam: File) extends CommandLineFunction with LongRunTime {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val xhmmParams = xhmmParamsArg

    @Input(doc = "")
    val origRD = origRDParam

    @Output
    val xcnv: File = new File(outputBase.getPath + ".xcnv")

    @Output
    val aux_xcnv: File = new File(outputBase.getPath + ".aux_xcnv")

    val posteriorsBase = outputBase.getPath

    @Output
    val dipPosteriors: File = new File(posteriorsBase + ".posteriors.DIP.txt")

    @Output
    val delPosteriors: File = new File(posteriorsBase + ".posteriors.DEL.txt")

    @Output
    val dupPosteriors: File = new File(posteriorsBase + ".posteriors.DUP.txt")

    var command: String =
      xhmmExec + " --discover" +
      " -p " + xhmmParams +
      " -r " + input +
      " -R " + origRD +
      " -c " + xcnv +
      " -a " + aux_xcnv +
      " -s " + posteriorsBase +
      " " + discoverCommandLineParams

    def commandLine = command

    override def description = "Discovers CNVs in normalized data: " + command
  }

  abstract class BaseGenotypeCNVs(inputParam: File, xcnv: File, origRDParam: File) extends CommandLineFunction with LongRunTime {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val xhmmParams = xhmmParamsArg

    @Input(doc = "")
    val origRD = origRDParam

    @Input(doc = "")
    val inXcnv = xcnv

    var command: String =
      xhmmExec + " --genotype" +
      " -p " + xhmmParams +
      " -r " + input +
      " -g " + inXcnv +
      " -F " + referenceFile +
      " -R " + origRD +
      " " + genotypeCommandLineParams
  }

  class GenotypeCNVs(inputParam: File, xcnv: File, origRDParam: File) extends BaseGenotypeCNVs(inputParam, xcnv, origRDParam) {
    @Output
    val vcf: File = new File(outputBase.getPath + ".vcf")

    command +=
      " -v " +  vcf

    def commandLine = command

    override def description = "Genotypes discovered CNVs in all samples: " + command
  }

  class GenotypeCNVandSubsegments(inputParam: File, xcnv: File, origRDParam: File) extends BaseGenotypeCNVs(inputParam, xcnv, origRDParam) {
    @Output
    val vcf: File = new File(outputBase.getPath + ".subsegments.vcf")

    command +=
      " -v " +  vcf +
      " --subsegments" +
      " --maxTargetsInSubsegment " + maxTargetsInSubsegment +
      " --genotypeQualThresholdWhenNoExact " + subsegmentGenotypeThreshold

    def commandLine = command

    override def description = "Genotypes discovered CNVs (and their sub-segments, of up to " + maxTargetsInSubsegment + " targets) in all samples: " + command
  }
}
