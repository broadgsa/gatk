/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.extensions.gatk

import java.io.{PrintStream, PrintWriter, File}
import org.broadinstitute.gatk.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.gatk.utils.commandline.{Input, Gather, Output}
import org.broadinstitute.gatk.queue.function.{InProcessFunction, CommandLineFunction}
import org.broadinstitute.gatk.tools.walkers.coverage.CoverageUtils
import scala.collection.JavaConversions._
import scala.Some
import org.broadinstitute.gatk.utils.text.XReadLines
import org.broadinstitute.gatk.queue.util.VCF_BAM_utilities
import org.broadinstitute.gatk.utils.downsampling.DownsampleType

// Minimal refactor from a package object to a file full of classes/objects
// due to ongoing bugs with inner classes/objects in package objects:
//   https://issues.scala-lang.org/browse/SI-4344
//   https://issues.scala-lang.org/browse/SI-5954

  class DoC(val bams: List[File], val DoC_output: File, val countType: CoverageUtils.CountPileupType, val MAX_DEPTH: Int, val minMappingQuality: Int, val minBaseQuality: Int, val scatterCountInput: Int, val START_BIN: Int, val NUM_BINS: Int, val minCoverageCalcs: Seq[Int], val sampleRenameMappingFile: Option[File] = None) extends CommandLineGATK with ScatterGatherableFunction {
    val DOC_OUTPUT_SUFFIX: String = ".sample_interval_summary"

    // So that the output files of this DoC run get deleted once they're used further downstream:
    this.isIntermediate = true

    this.analysis_type = "DepthOfCoverage"

    this.input_file = bams
    if (sampleRenameMappingFile.isDefined)
      this.sample_rename_mapping_file = sampleRenameMappingFile.get

    this.downsample_to_coverage = Some(MAX_DEPTH)
    this.downsampling_type = DownsampleType.BY_SAMPLE

    this.scatterCount = scatterCountInput
    this.scatterClass = classOf[IntervalScatterFunction]

    // HACK for DoC to work properly within Queue:
    @Output
    @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
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

  class DoCwithDepthOutputAtEachBase(bams: List[File], DoC_output: File, countType: CoverageUtils.CountPileupType, MAX_DEPTH: Int, minMappingQuality: Int, minBaseQuality: Int, scatterCountInput: Int, START_BIN: Int, NUM_BINS: Int, minCoverageCalcs: Seq[Int], sampleRenameMappingFile: Option[File] = None) extends DoC(bams, DoC_output, countType: CoverageUtils.CountPileupType, MAX_DEPTH: Int, minMappingQuality, minBaseQuality, scatterCountInput, START_BIN, NUM_BINS, minCoverageCalcs, sampleRenameMappingFile) {
    // HACK for DoC to work properly within Queue:
    @Output
    @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
    var outPrefix = DoC_output

    override def commandLine = super.commandLine.replaceAll(" --omitDepthOutputAtEachBase", "")
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

  class MergeGATKdepths(DoCsToCombine: List[File], outFile: String, columnSuffix: String, xhmmExec: File, rdPrecisionArg: Option[Int], outputTargetsBySamples: Boolean, sampleIDsMap: String = "", sampleIDsMapFromColumn: Int = 1, sampleIDsMapToColumn: Int = 2) extends CommandLineFunction {
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

  class ParsedBamListWithOptionalSampleMappings(bamsFile: File) {
    var bams = bamsFile

    var allBams = List[File]()
    var bamsWithoutSampleMapping = List[File]()
    var userMappedSampleToBams = scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]

    var sampleToBams = scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]
  }

object DoC {
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

  def parseBamListWithOptionalSampleMappings(bamsFile: File): ParsedBamListWithOptionalSampleMappings = {
    val r = new ParsedBamListWithOptionalSampleMappings(bamsFile)

    val rows = asScalaIterator(new XReadLines(r.bams))

    while (rows.hasNext) {
      val line = rows.next
      val splitLine = line.split("\\t")

      if (splitLine.length < 1 || splitLine.length > 2)
        throw new Exception("Invalid row in " + bamsFile.getPath + " : " + line)

      val bam = splitLine(0)
      val bamFile = new File(bam)
      r.allBams ::= bamFile

      if (splitLine.length == 2) {
        val sampleName = splitLine(1)

        if (r.userMappedSampleToBams.contains(sampleName))
          throw new Exception("Cannot map multiple BAM files to the same sample name: " + sampleName)

        r.userMappedSampleToBams += sampleName -> (scala.collection.mutable.Set.empty[File] + bamFile)
      }
      else {
        r.bamsWithoutSampleMapping ::= bamFile
      }
    }

    val autoMappedSampleToBams = VCF_BAM_utilities.getMapOfBAMsForSample(r.bamsWithoutSampleMapping)

    val overlappingSamples = autoMappedSampleToBams.keys.toList.intersect(r.userMappedSampleToBams.keys.toList)
    if (overlappingSamples.nonEmpty)
      throw new Exception("Cannot have the same sample mapped to different BAMs: " + overlappingSamples.toString)

    r.sampleToBams = autoMappedSampleToBams
    r.userMappedSampleToBams.foreach{ keyVal => {r.sampleToBams += keyVal._1 -> keyVal._2} }


    return r
  }
}

  class ProcessBamListWithOptionalSampleMappings(parsedBamList: ParsedBamListWithOptionalSampleMappings, outputBase: String) extends InProcessFunction {
    @Input(doc="")
    var bams: File = parsedBamList.bams

    @Output(doc="")
    var bamsList: File = new File(outputBase + ".bam.list")

    @Output(doc="")
    var bamSampleMap: File = new File(outputBase + ".bam_sample.txt")

    def run = {
      val bamsListWriter = new PrintWriter(new PrintStream(bamsList))
      for (bam <- parsedBamList.allBams) {
        bamsListWriter.printf("%s%n", bam)
      }
      bamsListWriter.close

      val bamSampleMapWriter = new PrintWriter(new PrintStream(bamSampleMap))
      for ((sampleName, sampleNameBams) <- parsedBamList.userMappedSampleToBams) {
        sampleNameBams.foreach { bam => bamSampleMapWriter.printf("%s\t%s%n", bam, sampleName) }
      }
      bamSampleMapWriter.close
    }
  }
