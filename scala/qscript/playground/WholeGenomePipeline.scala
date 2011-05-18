/*
 * Copyright (c) 2011, The Broad Institute 
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils

class WholeGenomePipeline extends QScript {
  @Input(doc="Bam file list", shortName = "I", required=true)
  var bamList: File = _

  @Input(doc="Exclude intervals list", shortName = "XL", required=false)
  var excludeIntervals: List[File] = Nil

  @Argument(doc="path to tmp space for storing intermediate bam files", shortName="cleanerTmpDir", required=true)
  var cleanerTmpDir: String = _

  @Argument(doc="Flag for running the whole genome (wg) or chromosome 20 (chr20). The default is chr20.", shortName="runType", required=false)
  var runType = "chr20"

  @Argument(doc="Chunk size. Defaults to 3,000,000", shortName="chunk", required=false)
  var chunkSize = 3000000

  val resources = "/humgen/gsa-pipeline/resources/5777/b37/"
  val reference = resources + "human_g1k_v37.fasta"
  val dbsnp = resources + "dbsnp_132.b37.vcf"
  val indels = resources + "1000G_indels_for_realignment.b37.vcf"
  val omni = resources + "1000G_omni2.5.b37.sites.vcf"
  val hapmap = resources + "hapmap_3.3.b37.sites.vcf"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.intervalsString = runIntervals
    this.memoryLimit = 2
  }

  case class ChrRange(chr: String, start: Long, stop: Long) {
    def toInterval = "%s:%d-%d".format(chr, start, stop)
    override def toString = toInterval
  }

  var runIntervals = List.empty[String]

  def script() {

    val rangeMap = Map(
      "cent1" -> new ChrRange("1", 121429168, 121529168),
      "cent16" -> new ChrRange("16", 40844464, 40944464),
      "chr20" -> new ChrRange("20", 1, 63025520),
      "chr20_100k" -> new ChrRange("20", 100001, 200000))

    var ranges = Traversable.empty[ChrRange]
    val chrs = IntervalUtils.getContigSizes(reference)

    runType = runType.toLowerCase
    if (runType == "wg") {
        // use all chromosomes from 1 to their length
        ranges = chrs.map { case (chr, len) => new ChrRange(chr, 1, len.longValue) }
        runIntervals = Nil
    } else {
      rangeMap.get(runType) match {
        case Some(range) =>
          ranges = List(range)
          runIntervals = List(range.toInterval)
        case None =>
          throw new RuntimeException("Inalid runType: " + runType + ". Must be one of: " + rangeMap.keys.mkString(", ") + ", or wg")
      }
    }

    val project = Array(".bams.list", ".bam.list", ".list").foldLeft(bamList.getName)(_.stripSuffix(_))
    val projectBase = project + "." + runType

    var chunkVcfs = List.empty[File]
    for (range <- ranges) {
      val chr = range.chr
      val chrStart = range.start
      val chrStop = range.stop

      var start = chrStart
      var chunkNumber = 1

      while (start <= chrStop) {
        val stop = (start + chunkSize - 1) min chrStop

        val chunkBase: String = "chunks/" + project + "." + runType + "_chunk_" + chr + "_" + chunkNumber
        val tmpBase: String = cleanerTmpDir + "/" + chunkBase

        val chunkInterval = List("%s:%d-%d".format(chr, start, stop))

        val target = new RealignerTargetCreator with CommandLineGATKArgs
        target.input_file :+= bamList
        target.intervalsString = chunkInterval
        target.excludeIntervals = excludeIntervals
        target.mismatchFraction = 0.0
        target.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        target.rodBind :+= RodBind("indels", "VCF", indels)
        target.out = tmpBase + ".target.intervals"
        target.jobOutputFile = chunkBase + ".target.intervals.out"
        target.isIntermediate = true
        add(target)

        val clean = new IndelRealigner with CommandLineGATKArgs
        clean.input_file :+= bamList
        clean.intervalsString = chunkInterval
        clean.excludeIntervals = excludeIntervals
        clean.targetIntervals = target.out
        clean.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        clean.rodBind :+= RodBind("indels", "VCF", indels)
        clean.doNotUseSW = true
        clean.simplifyBAM = true
        clean.bam_compression = 1
        clean.out = tmpBase + ".cleaned.bam"
        clean.jobOutputFile = chunkBase + ".cleaned.bam.out"
        clean.isIntermediate = true
        add(clean)

        val call = new UnifiedGenotyper with CommandLineGATKArgs
        call.input_file :+= clean.out
        call.intervalsString = chunkInterval
        call.excludeIntervals = excludeIntervals
        call.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        call.downsample_to_coverage = 50
        call.standard_min_confidence_threshold_for_calling = 4.0
        call.standard_min_confidence_threshold_for_emitting = 4.0
        call.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
        call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
        call.GSA_PRODUCTION_ONLY = true
        call.out = chunkBase + ".vcf"
        call.jobOutputFile = call.out + ".out"
        add(call)

        chunkVcfs :+= call.out
        start += chunkSize
        chunkNumber += 1
      }
    }

    val combineVCFs = new CombineVariants with CommandLineGATKArgs
    combineVCFs.rodBind = chunkVcfs.zipWithIndex.map { case (vcf, index) => RodBind("input"+index, "VCF", vcf) }
    combineVCFs.rod_priority_list = chunkVcfs.zipWithIndex.map { case (vcf, index) => "input"+index }.mkString(",")
    combineVCFs.variantmergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION
    combineVCFs.genotypemergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.PRIORITIZE
    combineVCFs.out = projectBase + ".merged.vcf"
    combineVCFs.jobOutputFile = combineVCFs.out + ".out"
    combineVCFs.assumeIdenticalSamples = true
    add(combineVCFs)

    val sv = new SelectVariants with CommandLineGATKArgs
    sv.selectIndels = true
    sv.rodBind :+= RodBind("variant", "VCF", combineVCFs.out)
    sv.out = projectBase + ".indels.vcf"
    sv.jobOutputFile = sv.out + ".out"
    add(sv)

    val filter = new VariantFiltration with CommandLineGATKArgs
    filter.variantVCF = sv.out
    filter.filterName :+= "HARD_TO_VALIDATE"
    filter.filterExpression :+= "\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\""
    filter.out = swapExt(sv.out, ".vcf", ".filtered.vcf")
    filter.jobOutputFile = filter.out + ".out"
    add(filter)

    val recombine = new CombineVariants with CommandLineGATKArgs
    recombine.rodBind :+= RodBind("indels", "VCF", sv.out)
    recombine.rodBind :+= RodBind("all", "VCF", combineVCFs.out)
    recombine.rod_priority_list = "indels,all"
    recombine.genotypemergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.PRIORITIZE
    recombine.out = swapExt(combineVCFs.out, ".vcf", ".filtered.vcf")
    recombine.jobOutputFile = recombine.out + ".out"
    add(recombine)

    val vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.rodBind :+= RodBind("input", "VCF", recombine.out)
    vr.rodBind :+= RodBind("hapmap", "VCF", hapmap, "known=false,training=true,truth=true,prior=15.0")
    vr.rodBind :+= RodBind("omni", "VCF", omni, "known=false,training=true,truth=false,prior=12.0")
    vr.rodBind :+= RodBind("dbsnp", "VCF", dbsnp, "known=true,training=false,truth=false,prior=10.0")
    vr.use_annotation = List("QD", "SB", "HaplotypeScore", "HRun")
    vr.trustAllPolymorphic = true
    vr.TStranche = List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    vr.tranches_file = swapExt(filter.out, ".vcf", ".tranches")
    vr.recal_file = swapExt(filter.out, ".vcf", ".recal")
    vr.jobOutputFile = vr.recal_file + ".out"
    add(vr)

    val ar = new ApplyRecalibration with CommandLineGATKArgs
    ar.rodBind :+= RodBind("input", "VCF", recombine.out)
    ar.tranches_file = vr.tranches_file
    ar.recal_file = vr.recal_file
    ar.ts_filter_level = 99.0
    ar.ignore_filter :+= "HARD_TO_VALIDATE"
    ar.out = swapExt(recombine.out, ".vcf", ".recalibrated.vcf")
    ar.jobOutputFile = ar.out + ".out"
    add(ar)

    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.tranchesFile = vr.tranches_file
    stdEval.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
    stdEval.rodBind :+= RodBind("eval", "VCF", ar.out)
    stdEval.doNotUseAllStandardStratifications = true
    stdEval.doNotUseAllStandardModules = true
    stdEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    stdEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
    stdEval.out = swapExt(ar.out, ".vcf", ".eval")
    stdEval.jobOutputFile = stdEval.out + ".out"
    add(stdEval)
  }
}
