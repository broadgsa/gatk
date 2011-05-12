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

  @Argument(doc="path to tmp space for storing intermediate bam files", shortName="cleanerTmpDir", required=true)
  var cleanerTmpDir: String = _

  @Argument(doc="Flag for running the whole genome (wg), chromosome 20 (chr20) or just the centromere of chromosome 1 (cent1). The default is cent1.", shortName="runType", required=false)
  var runType = "cent1"

  @Argument(doc="Chunk size. Defaults to 3,000,000", shortName="chunk", required=false)
  var chunkSize = 3000000

  val reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
  val dbsnp = "/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf"
  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.intervalsString = runIntervals
    this.memoryLimit = 2
  }

  case class ChrRange(chr: String, start: Long, stop: Long) {
    def toInterval = "%s:%d-%d".format(chr, start, stop)
    override def toString = toInterval
  }

  val chr1Centromere = new ChrRange("1", 121429168, 121529169)

  var runIntervals = List.empty[String]

  def script() {

    var ranges = Traversable.empty[ChrRange]
    val chrs = IntervalUtils.getContigSizes(reference)
    runType match {
      case "wg" =>
        // use all chromosomes from 1 to their length
        ranges = chrs.map { case (chr, len) => new ChrRange(chr, 1, len.longValue) }
        runIntervals = Nil
      case "chr20" =>
        // use chromosome 20 from 1 to its length
        ranges = chrs.filterKeys(_ == "20").map { case (chr, len) => new ChrRange(chr, 1, len.longValue) }
        runIntervals = List(ranges.head.toInterval)
      case "cent1" =>
        // use chromosome 1 centromere
        ranges = List(chr1Centromere)
        runIntervals = List(ranges.head.toInterval)
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
        var interval = "%s:%d-%d".format(chr, start, stop)
        val chunkBase: String = "chunks/" + project + "." + runType + "_chunk_" + chr + "_" + chunkNumber
        val tmpBase: String = cleanerTmpDir + "/" + chunkBase

        val target = new RealignerTargetCreator with CommandLineGATKArgs
        target.intervalsString :+= interval
        target.input_file :+= bamList
        target.mismatchFraction = 0.0
        target.maxIntervalSize = 700
        target.out = tmpBase + ".target.intervals"
        target.jobOutputFile = chunkBase + ".target.out"
        target.isIntermediate = true
        add(target)

        val clean = new IndelRealigner with CommandLineGATKArgs
        clean.intervalsString :+= interval
        clean.input_file :+= bamList
        clean.targetIntervals = target.out
        clean.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        clean.doNotUseSW = true
        clean.simplifyBAM = true
        clean.bam_compression = 1
        clean.out = tmpBase + ".cleaned.bam"
        clean.jobOutputFile = chunkBase + ".clean.out"
        clean.jobPriority = 51
        clean.isIntermediate = true
        add(clean)

        val call = new UnifiedGenotyper with CommandLineGATKArgs
        call.intervalsString :+= interval
        call.input_file :+= clean.out
        call.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        call.downsample_to_coverage = 50
        call.standard_min_confidence_threshold_for_calling = 4.0
        call.standard_min_confidence_threshold_for_emitting = 4.0
        call.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
        call.exactCalculation = org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel.ExactCalculation.LINEAR_EXPERIMENTAL
        call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
        call.GSA_PRODUCTION_ONLY = true
        call.out = chunkBase + ".vcf"
        call.jobOutputFile = chunkBase + ".call.out"
        call.jobPriority = 52
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
    combineVCFs.jobOutputFile = projectBase + ".combine.out"
    combineVCFs.assumeIdenticalSamples = true
    add(combineVCFs)

    val sv = new SelectVariants with CommandLineGATKArgs
    sv.selectIndels = true
    sv.rodBind :+= RodBind("variant", "VCF", combineVCFs.out)
    sv.out = swapExt(combineVCFs.out, ".merged.vcf", ".indels.vcf")
    sv.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".sv.out")
    add(sv)

    val filter = new VariantFiltration with CommandLineGATKArgs
    filter.variantVCF = sv.out
    filter.filterName :+= "HARD_TO_VALIDATE"
    filter.filterExpression :+= "\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\""
    filter.out = swapExt(sv.out, ".vcf", ".filtered.vcf")
    filter.jobOutputFile = swapExt(sv.jobOutputFile, ".sv.out", ".filter.out")
    add(filter)

    val recombine = new CombineVariants with CommandLineGATKArgs
    recombine.rodBind :+= RodBind("indels", "VCF", sv.out)
    recombine.rodBind :+= RodBind("all", "VCF", combineVCFs.out)
    recombine.rod_priority_list = "indels,all"
    recombine.genotypemergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.PRIORITIZE
    recombine.out = swapExt(combineVCFs.out, ".vcf", ".filtered.vcf")
    recombine.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".recombine.out")
    add(recombine)

    val vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.rodBind :+= RodBind("input", "VCF", recombine.out)
    vr.rodBind :+= RodBind("hapmap", "VCF", hapmap, "training=true, prior=3.0")
    vr.rodBind :+= RodBind("1kg", "VCF", omni, "training=true, prior=3.0")
    vr.rodBind :+= RodBind("dbsnp", "VCF", dbsnp, "known=true")
    vr.use_annotation = List("QD", "SB", "HaplotypeScore", "HRun")
    vr.trustAllPolymorphic = true
    vr.TStranche = List("0.1","0.5","0.75","1.0","1.25","1.35","1.5","1.7","1.8","1.9","2.0","2.1","2.2","2.3","2.5","3.0","5.0","8.0","10.0")
    vr.tranches_file = swapExt(filter.out, ".merged.filtered.vcf", ".tranches")
    vr.recal_file = swapExt(filter.out, ".merged.filtered.vcf", "")
    vr.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".cr.out")
    add(vr)

    val ar = new ApplyRecalibration with CommandLineGATKArgs
    ar.tranches_file = vr.tranches_file
    ar.recal_file = vr.recal_file
    ar.ignore_filter :+= "HARD_TO_VALIDATE"
    ar.out = swapExt(filter.out, ".vcf", ".recalibrated.vcf")
    ar.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".ar.out")
    add(ar)

    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.tranchesFile = vr.tranches_file
    stdEval.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
    stdEval.rodBind :+= RodBind("eval", "VCF", ar.out)
    stdEval.doNotUseAllStandardStratifications = true
    stdEval.doNotUseAllStandardModules = true
    stdEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    stdEval.stratificationModule = List("EvalRod", "CompRod", "Novelty")
    stdEval.out = swapExt(ar.out, ".vcf", ".eval")
    stdEval.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".eval.out")
    add(stdEval)
  }
}
