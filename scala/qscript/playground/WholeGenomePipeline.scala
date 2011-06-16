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

import io.Source
import org.broadinstitute.sting.queue.extensions.samtools.{SamtoolsIndexFunction, SamtoolsMergeFunction}
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
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

  @Argument(doc="Memory limit. Defaults to 4g", shortName="pipeMemory", required=false)
  var pipelineMemoryLimit = 4

  val resources = "/humgen/gsa-pipeline/resources/5777/b37/"
  val reference = resources + "human_g1k_v37.fasta"
  val dbsnp = resources + "dbsnp_132.b37.vcf"
  val indels = resources + "1000G_indels_for_realignment.b37.vcf"
  val omni = resources + "1000G_omni2.5.b37.sites.vcf"
  val hapmap = resources + "hapmap_3.3.b37.sites.vcf"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.intervalsString = runIntervals
    this.memoryLimit = pipelineMemoryLimit
  }

  case class Interval(chr: String, start: Long, stop: Long) {
    override def toString = "%s:%d-%d".format(chr, start, stop)
  }

  var runIntervals = List.empty[String]

  def script() {

    var intervals = Traversable.empty[Interval]

    runType = runType.toLowerCase
    if (runType == "wg") {
      val contigs = (1 to 22).map(_.toString) ++ List("X", "Y", "MT")
      val sizes = IntervalUtils.getContigSizes(reference)
      intervals = contigs.map(chr => new Interval(chr, 1, sizes.get(chr).longValue))
      runIntervals = Nil
    } else {
      val locs = Map(
        "cent1" -> new Interval("1", 121429168, 121529168),
        "cent16" -> new Interval("16", 40844464, 40944464),
        "chr20" -> new Interval("20", 1, 63025520),
        "chr20_100k" -> new Interval("20", 100001, 200000))

      locs.get(runType) match {
        case Some(range) =>
          intervals = List(range)
          runIntervals = List(range.toString)
        case None =>
          throw new RuntimeException("Invalid runType: " + runType + ". Must be one of: " + locs.keys.mkString(", ") + ", or wg")
      }
    }

    val project = Array(".bams.list", ".bam.list", ".list").foldLeft(bamList.getName)(_.stripSuffix(_))
    val projectBase = project + "." + runType

    val mergeBam = new SamtoolsMergeFunction
    mergeBam.inputBams = Source.fromFile(bamList).getLines().toList
    if (runType != "wg")
      mergeBam.region = intervals.head.toString
    mergeBam.memoryLimit = pipelineMemoryLimit
    mergeBam.outputBam = cleanerTmpDir + "/" + projectBase + ".unclean.bam"
    mergeBam.jobOutputFile = projectBase + ".unclean.bam.out"
    mergeBam.isIntermediate = true
    mergeBam.memoryLimit = pipelineMemoryLimit
    add(mergeBam)

    val indexBam = new SamtoolsIndexFunction
    indexBam.bamFile = mergeBam.outputBam
    indexBam.memoryLimit = pipelineMemoryLimit
    indexBam.jobOutputFile = projectBase + ".unclean.bam.bai.out"
    indexBam.isIntermediate = true
    add(indexBam)

    var chunkVcfs = List.empty[File]
    for (interval <- intervals) {
      val chr = interval.chr
      val chrStart = interval.start
      val chrStop = interval.stop

      var start = chrStart
      var chunkNumber = 1

      while (start <= chrStop) {
        val stop = (start + chunkSize - 1) min chrStop

        val chunkBase: String = "chunks/" + project + "." + runType + "_chunk_" + chr + "_" + chunkNumber
        val tmpBase: String = cleanerTmpDir + "/" + chunkBase

        val chunkInterval = List("%s:%d-%d".format(chr, start, stop))

        val target = new RealignerTargetCreator with CommandLineGATKArgs
        target.input_file :+= mergeBam.outputBam
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
        clean.input_file :+= mergeBam.outputBam
        clean.intervalsString = chunkInterval
        clean.excludeIntervals = excludeIntervals
        clean.targetIntervals = target.out
        clean.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
        clean.rodBind :+= RodBind("indels", "VCF", indels)
        clean.consensusDeterminationModel = org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
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

    val combineChunks = new CombineVariants with CommandLineGATKArgs
    combineChunks.rodBind = chunkVcfs.zipWithIndex.map { case (vcf, index) => RodBind("input"+index, "VCF", vcf) }
    combineChunks.rod_priority_list = chunkVcfs.zipWithIndex.map { case (vcf, index) => "input"+index }.mkString(",")
    combineChunks.filteredrecordsmergetype = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineChunks.assumeIdenticalSamples = true
    combineChunks.out = projectBase + ".unfiltered.vcf"
    combineChunks.jobOutputFile = combineChunks.out + ".out"
    add(combineChunks)

    val selectSNPs = new SelectVariants with CommandLineGATKArgs
    selectSNPs.selectSNPs = true
    selectSNPs.rodBind :+= RodBind("variant", "VCF", combineChunks.out)
    selectSNPs.out = projectBase + ".snps.unrecalibrated.vcf"
    selectSNPs.jobOutputFile = selectSNPs.out + ".out"
    add(selectSNPs)

    val selectIndels = new SelectVariants with CommandLineGATKArgs
    selectIndels.selectIndels = true
    selectIndels.rodBind :+= RodBind("variant", "VCF", combineChunks.out)
    selectIndels.out = projectBase + ".indels.unfiltered.vcf"
    selectIndels.jobOutputFile = selectIndels.out + ".out"
    add(selectIndels)

    val filterIndels = new VariantFiltration with CommandLineGATKArgs
    filterIndels.variantVCF = selectIndels.out
    filterIndels.filterName = List("Indel_QUAL", "Indel_SB", "Indel_QD", "Indel_HRun", "Indel_HaplotypeScore")
    filterIndels.filterExpression = List("\"QUAL<30.0\"", "\"SB>-1.0\"", "\"QD<2.0\"", "\"HRun>15\"", "\"HaplotypeScore>20.0\"")
    filterIndels.out = projectBase + ".indels.filtered.vcf"
    filterIndels.jobOutputFile = filterIndels.out + ".out"
    add(filterIndels)

    val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs
    combineSNPsIndels.rodBind :+= RodBind("indels", "VCF", selectIndels.out)
    combineSNPsIndels.rodBind :+= RodBind("snps", "VCF", selectSNPs.out)
    combineSNPsIndels.rod_priority_list = "indels,snps"
    combineSNPsIndels.filteredRecordsMergeType = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineSNPsIndels.assumeIdenticalSamples = true
    combineSNPsIndels.out = projectBase + ".unrecalibrated.vcf"
    combineSNPsIndels.jobOutputFile = combineSNPsIndels.out + ".out"
    add(combineSNPsIndels)
 
    val vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.rodBind :+= RodBind("input", "VCF", combineSNPsIndels.out)
    vr.rodBind :+= RodBind("hapmap", "VCF", hapmap, "known=false,training=true,truth=true,prior=15.0")
    vr.rodBind :+= RodBind("omni", "VCF", omni, "known=false,training=true,truth=false,prior=12.0")
    vr.rodBind :+= RodBind("dbsnp", "VCF", dbsnp, "known=true,training=false,truth=false,prior=8.0")
    vr.trustAllPolymorphic = true
    vr.use_annotation = List("QD", "HaplotypeScore", "HRun", "MQRankSum", "ReadPosRankSum")
    vr.TStranche = List(
      "100.0", "99.9", "99.5", "99.3",
      "99.0", "98.9", "98.8",
      "98.5", "98.4", "98.3", "98.2", "98.1",
      "98.0",
      "97.0",
      "95.0",
      "90.0")
    vr.tranches_file = projectBase + ".tranches"
    vr.recal_file = projectBase + ".recal"
    vr.jobOutputFile = vr.recal_file + ".out"
    vr.memoryLimit = 32
    add(vr)

    for (tranche <- vr.TStranche) {
      val ar = new ApplyRecalibration with CommandLineGATKArgs
      ar.rodBind :+= RodBind("input", "VCF", combineSNPsIndels.out)
      ar.tranches_file = vr.tranches_file
      ar.recal_file = vr.recal_file
      ar.ts_filter_level = tranche.toDouble
      ar.out = projectBase + ".recalibrated." + tranche + ".vcf"
      ar.jobOutputFile = ar.out + ".out"
      ar.memoryLimit = 32
      add(ar)

      val eval = new VariantEval with CommandLineGATKArgs
      eval.tranchesFile = vr.tranches_file
      eval.rodBind :+= RodBind("eval", "VCF", ar.out)
      eval.rodBind :+= RodBind("dbsnp", "VCF", dbsnp)
      eval.doNotUseAllStandardStratifications = true
      eval.doNotUseAllStandardModules = true
      eval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
      eval.stratificationModule = List("EvalRod", "CompRod", "Novelty")
      eval.out = swapExt(ar.out, ".vcf", ".eval")
      eval.jobOutputFile = eval.out + ".out"
      eval.memoryLimit = 32
      add(eval)
    }
  }
}
