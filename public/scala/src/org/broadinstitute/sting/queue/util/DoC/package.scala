package org.broadinstitute.sting.queue.util

import java.io.File
import org.broadinstitute.sting.queue.extensions.gatk.{IntervalScatterFunction, CommandLineGATK}
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.gatk.downsampling.DownsampleType
import org.broadinstitute.sting.commandline.{Input, Gather, Output}
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.gatk.walkers.coverage.CoverageUtils

package object DoC {
  class DoC(val bams: List[File], val DoC_output: File, val countType: CoverageUtils.CountPileupType, val MAX_DEPTH: Int, val minMappingQuality: Int, val minBaseQuality: Int, val scatterCountInput: Int, val START_BIN: Int, val NUM_BINS: Int, val minCoverageCalcs: Seq[Int]) extends CommandLineGATK with ScatterGatherableFunction {
    val DOC_OUTPUT_SUFFIX: String = ".sample_interval_summary"

    // So that the output files of this DoC run get deleted once they're used further downstream:
    this.isIntermediate = true

    this.analysis_type = "DepthOfCoverage"

    this.input_file = bams

    this.downsample_to_coverage = Some(MAX_DEPTH)
    this.downsampling_type = DownsampleType.BY_SAMPLE

    this.scatterCount = scatterCountInput
    this.scatterClass = classOf[IntervalScatterFunction]

    // HACK for DoC to work properly within Queue:
    @Output
    @Gather(classOf[org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction])
    var intervalSampleOut: File = new File(DoC_output.getPath + DOC_OUTPUT_SUFFIX)

    override def commandLine = super.commandLine +
      " --omitDepthOutputAtEachBase" +
      " --omitLocusTable" +
      " --minMappingQuality " + minMappingQuality +
      " --minBaseQuality " + minBaseQuality +
      optional("--countType", countType, spaceSeparated=true, escape=true, format="%s") +
      " --start " + START_BIN + " --stop " + MAX_DEPTH + " --nBins " + NUM_BINS +
      (if (!minCoverageCalcs.isEmpty) minCoverageCalcs.map(cov => " --summaryCoverageThreshold " + cov).reduceLeft(_ + "" + _) else "") +
      " --includeRefNSites" +
      " -o " + DoC_output

    override def shortDescription = "DoC: " + DoC_output
  }

  class DoCwithDepthOutputAtEachBase(bams: List[File], DoC_output: File, countType: CoverageUtils.CountPileupType, MAX_DEPTH: Int, minMappingQuality: Int, minBaseQuality: Int, scatterCountInput: Int, START_BIN: Int, NUM_BINS: Int, minCoverageCalcs: Seq[Int]) extends DoC(bams, DoC_output, countType: CoverageUtils.CountPileupType, MAX_DEPTH: Int, minMappingQuality, minBaseQuality, scatterCountInput, START_BIN, NUM_BINS, minCoverageCalcs) {
    // HACK for DoC to work properly within Queue:
    @Output
    @Gather(classOf[org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction])
    var outPrefix = DoC_output

    override def commandLine = super.commandLine.replaceAll(" --omitDepthOutputAtEachBase", "")
  }

  def buildDoCgroups(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]], samplesPerJob: Int, outputBase: File): List[Group] = {
    var l: List[Group] = Nil

    var remaining = samples
    var subsamples: List[String] = Nil
    var count = 1

    while (!remaining.isEmpty) {
      val splitRes = (remaining splitAt samplesPerJob)
      subsamples = splitRes._1
      remaining = splitRes._2
      l ::= new Group("group" + count, outputBase, subsamples, VCF_BAM_utilities.findBAMsForSamples(subsamples, sampleToBams))
      count = count + 1
    }

    return l
  }

  // A group has a list of samples and bam files to use for DoC
  class Group(val name: String, val outputBase: File, val samples: List[String], val bams: List[File]) {
    // getName() just includes the file name WITHOUT the path:
    val groupOutputName = name + "." + outputBase.getName

    // Comment this out, so that when jobs are scattered in DoC class below, they do not scatter into outputs at directories that don't exist!!! :
    //def DoC_output = new File(outputBase.getParentFile(), groupOutputName)

    def DoC_output = new File(groupOutputName)

    override def toString(): String = String.format("[Group %s [%s] with samples %s against bams %s]", name, DoC_output, samples, bams)
  }

  class MergeGATKdepths(DoCsToCombine: List[File], outFile: String, columnSuffix: String, xhmmExec: File, sampleIDsMap: String, sampleIDsMapFromColumn: Int, sampleIDsMapToColumn: Int, rdPrecisionArg: Option[Int], outputTargetsBySamples: Boolean) extends CommandLineFunction {
    @Input(doc = "")
    var inputDoCfiles: List[File] = DoCsToCombine

    @Output
    val mergedDoC: File = new File(outFile)
    var command: String =
      xhmmExec + " --mergeGATKdepths" +
        inputDoCfiles.map(input => " --GATKdepths " + input).reduceLeft(_ + "" + _) +
        " --columnSuffix " + columnSuffix +
        " -o " + mergedDoC
    if (sampleIDsMap != "")
      command += " --sampleIDmap " + sampleIDsMap + " --fromID " + sampleIDsMapFromColumn + " --toID " + sampleIDsMapToColumn
    rdPrecisionArg match {
      case Some(rdPrecision) => {
        command += " --rdPrecision " + rdPrecision
      }
      case None => {}
    }
    if (outputTargetsBySamples)
      command += " --outputTargetsBySamples"

    def commandLine = command

    override def description = "Combines DoC outputs for multiple samples (at same loci): " + command
  }

  class PrepareTargets(intervalsIn: List[File], outIntervals: String, val xhmmExec: File, val referenceFile: File) extends CommandLineFunction {
    @Input(doc = "List of files containing targeted intervals to be prepared and merged")
    var inIntervals: List[File] = intervalsIn

    @Output(doc = "The merged intervals file to write to")
    var out: File = new File(outIntervals)

    var command: String =
      xhmmExec + " --prepareTargets" +
        " -F " + referenceFile +
        inIntervals.map(intervFile => " --targets " + intervFile).reduceLeft(_ + "" + _) +
        " --mergedTargets " + out

    def commandLine = command

    override def description = "Sort all target intervals, merge overlapping ones, and print the resulting interval list: " + command
  }
}
