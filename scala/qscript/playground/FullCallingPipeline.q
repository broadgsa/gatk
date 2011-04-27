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

import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.broad.PicardPipeline

class FullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="level of parallelism for UnifiedGenotyper (both for SNPs and indels). By default is set to 20.", shortName="varScatter", required=false)
  var variantCallerScatterCount = 20

  @Argument(doc="memory limit for UnifiedGenotyper (both for SNPs and indels). By default is set to 4g.", shortName="varMemory", required=false)
  var variantCallerMemory = 4

  @Argument(doc="expand each target in input intervals by the specified number of bases (50 bases by default)", shortName="expand", required=false)
  var expandIntervals = 50

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.intervals = List(qscript.pipeline.getProject.getIntervalList)
    this.memoryLimit = 4
  }

  def script() {
    pipeline = PicardPipeline.parse(qscript.yamlFile)

    val projectBase: String = qscript.pipeline.getProject.getName
    val base = projectBase + ".cleaned"
    val bamType = "cleaned"

    // Make the bam list
    val writeBamList = new ListWriterFunction
    writeBamList.analysisName = base + "_BamList"
    writeBamList.jobOutputFile = ".queue/logs/Overall/WriteBamList.out"
    writeBamList.inputFiles = qscript.pipeline.getSamples.filter(_.getBamFiles.contains(bamType)).map(_.getBamFiles.get(bamType)).toList
    writeBamList.listFile = "Resources/" + base +".bamfiles.list"
    add(writeBamList)

    val ei = new ExpandIntervals(
      qscript.pipeline.getProject.getIntervalList,
      1,
      qscript.expandIntervals,
      "Resources/" + base + ".flanks.interval_list",
      qscript.pipeline.getProject.getReferenceFile,
      "INTERVALS",
      "INTERVALS")
    ei.jobOutputFile = ".queue/logs/Overall/ExpandIntervals.out"

    if (qscript.expandIntervals > 0) {
      add(ei)
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0) {
        this.intervals :+= ei.outList
      }
    }

    // Call indels
    val indels = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    indels.analysisName = base + "_indels"
    indels.jobOutputFile = ".queue/logs/IndelCalling/UnifiedGenotyper.indels.out"
    indels.downsample_to_coverage = 600
    indels.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    indels.input_file = List(writeBamList.listFile)
    indels.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
    indels.out = "IndelCalls/" + base+".indels.vcf"
    indels.scatterCount = qscript.variantCallerScatterCount
    indels.memoryLimit = qscript.variantCallerMemory
    add(indels)

    // Filter indels
    val filteredIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filteredIndels.analysisName = base + "_filteredIndels"
    filteredIndels.jobOutputFile = ".queue/logs/IndelCalling/VariantFiltration.indels.out"
    filteredIndels.filterName ++= List("IndelQUALFilter","IndelSBFilter","IndelQDFilter")
    filteredIndels.filterExpression ++= List("\"QUAL<30.0\"","\"SB>-1.0\"","\"QD<2.0\"")
    filteredIndels.variantVCF = indels.out
    filteredIndels.out = swapExt("IndelCalls", indels.out, ".vcf",".filtered.vcf")
    add(filteredIndels)

    // Call snps
    val snps = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    snps.analysisName = base+"_snps"
    snps.jobOutputFile = ".queue/logs/SNPCalling/UnifiedGenotyper.snps.out"
    snps.downsample_to_coverage = 600
    snps.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    snps.input_file = List(writeBamList.listFile)
    snps.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
    snps.out = "SnpCalls/" + base+".snps.vcf"
    snps.scatterCount = qscript.variantCallerScatterCount
    snps.memoryLimit = qscript.variantCallerMemory
    add(snps)

    // Filter snps
    val filteredSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filteredSNPs.analysisName = base+"_filteredSNPs"
    filteredSNPs.jobOutputFile = ".queue/logs/SNPCalling/VariantFiltration.snps.out"
    filteredSNPs.filterName ++= List("SNPSBFilter","SNPQDFilter","SNPHRunFilter")
    filteredSNPs.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    filteredSNPs.clusterWindowSize = 10
    filteredSNPs.clusterSize = 3
    filteredSNPs.rodBind :+= RodBind("mask", "VCF", filteredIndels.out)
    filteredSNPs.variantVCF = snps.out
    filteredSNPs.out = swapExt("SnpCalls",snps.out,".vcf",".filtered.vcf")
    add(filteredSNPs)

    // Combine indels and snps into one VCF
    val combineAll = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineAll.analysisName = base + "_combineAll"
    combineAll.jobOutputFile = ".queue/logs/Combined/CombineVariants.out"
    combineAll.variantmergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION
    combineAll.rod_priority_list = "Indels,SNPs"
    combineAll.rodBind :+= RodBind("Indels", "VCF", filteredIndels.out)
    combineAll.rodBind :+= RodBind("SNPs", "VCF", filteredSNPs.out)
    combineAll.out = "CombinedCalls/" + base + ".snps_and_indels.filtered.vcf"
    add(combineAll)

    // Annotate variants
    val annotated = new GenomicAnnotator with CommandLineGATKArgs with ExpandedIntervals
    annotated.analysisName = base+"_annotated"
    annotated.jobOutputFile = ".queue/logs/Combined/GenomicAnnotator.out"
    annotated.rodToIntervalTrackName = "variant"
    annotated.rodBind :+= RodBind("variant", "VCF", combineAll.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.pipeline.getProject.getRefseqTable)
    annotated.out = base + ".snps_and_indels.filtered.annotated.vcf"
    add(annotated)

    // Variant eval the standard region
    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.analysisName = base+"_StandardVariantEval"
    stdEval.jobOutputFile = ".queue/logs/Overall/VariantEval.std.out"
    stdEval.doNotUseAllStandardStratifications = true
    stdEval.doNotUseAllStandardModules = true
    stdEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    stdEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
    stdEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
    stdEval.rodBind :+= RodBind("eval", "VCF", annotated.out)
    stdEval.out = swapExt(annotated.out, ".vcf", ".eval")
    add(stdEval)

    if (qscript.expandIntervals > 0) {
      // Variant eval the flanking region
      val flanksEval = new VariantEval with CommandLineGATKArgs
      flanksEval.analysisName = base+"_FlanksVariantEval"
      flanksEval.jobOutputFile = ".queue/logs/Overall/VariantEval.flanks.out"
      flanksEval.intervals = List(ei.outList)
      flanksEval.doNotUseAllStandardStratifications = true
      flanksEval.doNotUseAllStandardModules = true
      flanksEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
      flanksEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
      flanksEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
      flanksEval.rodBind :+= RodBind("eval", "VCF", annotated.out)
      flanksEval.out = swapExt(annotated.out, ".vcf", ".flanks.eval")
      add(flanksEval)
    }
  }
}
