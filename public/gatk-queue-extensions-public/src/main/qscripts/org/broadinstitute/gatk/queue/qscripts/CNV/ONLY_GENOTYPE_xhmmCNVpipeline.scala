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

package org.broadinstitute.gatk.queue.qscripts.CNV

import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.util.VCF_BAM_utilities
import org.broadinstitute.gatk.queue.extensions.gatk.DoC._
import org.broadinstitute.gatk.utils.commandline._
import java.io.{File, PrintStream, PrintWriter}
import org.broadinstitute.gatk.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.gatk.tools.walkers.coverage.CoverageUtils
import org.broadinstitute.gatk.queue.function.scattergather.{CloneFunction, ScatterFunction, GatherFunction, ScatterGatherableFunction}
import org.broadinstitute.gatk.queue.function.{CommandLineFunction, InProcessFunction}
import org.broadinstitute.gatk.utils.io.IOUtils

class ONLY_GENOTYPE_xhmmCNVpipeline extends QScript {
  qscript =>

  @Input(doc = "bam input, as as a list of .bam files, or a list of bam files with sample IDs to be used ( as specified at https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--sample_rename_mapping_file )", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "xhmm executable file", shortName = "xhmmExec", required = true)
  var xhmmExec: File = _

  @Input(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(doc = "Samples to run together for DoC, CNV discovery, and CNV genotyping.  By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Base name for files to output", shortName = "o", required = true)
  var outputBase: File = _

  @Argument(shortName = "xhmmParams", doc = "xhmm model parameters file", required = true)
  var xhmmParamsArg: File = _

  @Argument(shortName = "genotypeParams", doc = "xhmm command-line parameters for genotyping step", required = false)
  var genotypeCommandLineParams: String = ""

  @Argument(shortName = "addGenotypeRegions", doc = "Additional interval list files to be genotyped", required = false)
  var addGenotypeRegions: List[File] = List[File]()

  @Argument(shortName = "longJobQueue", doc = "Job queue to run the 'long-running' commands", required = false)
  var longJobQueue: String = ""

  @Argument(shortName = "filteredZscored", doc = "File of PCA-normalized read depths, after filtering and Z-score calculation", required = true)
  var filteredZscored: File = _

  @Argument(shortName = "originalSameFiltered", doc = "File of original read depths, using same filters (samples and targets) as Z-score matrix [filteredZscored argument]", required = true)
  var originalSameFiltered: File = _


  trait LongRunTime extends CommandLineFunction {
    if (longJobQueue != "")
      this.jobQueue = longJobQueue
  }

  def script = {
    val parseMixedInputBamList = parseBamListWithOptionalSampleMappings(bams)

    val processMixedInputBamList = new ProcessBamListWithOptionalSampleMappings(parseMixedInputBamList, outputBase.getPath)
    add(processMixedInputBamList)

    val samples: List[String] = parseMixedInputBamList.sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val groups: List[Group] = buildDoCgroups(samples, parseMixedInputBamList.sampleToBams, samplesPerJob, outputBase)
    for (group <- groups) {
      Console.out.printf("Group is %s%n", group)
    }

    for (regionsFile <- addGenotypeRegions) {
    	val genotypeRegions = new GenotypeCNVs(filteredZscored, regionsFile, originalSameFiltered, new File(outputBase.getParent + "/" + regionsFile.getName.replace(".interval_list", "") + "." + outputBase.getName), xhmmParamsArg, referenceFile, genotypeCommandLineParams, xhmmExec, groups) with LongRunTime
	add(genotypeRegions)
    }
  }

}
