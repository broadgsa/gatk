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

class HybridSelectionPipeline extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="level of parallelism for UnifiedGenotyper. By default set to 20.", shortName="varScatter", required=false)
  var variantCallerScatterCount = 20

  @Argument(doc="memory limit for UnifiedGenotyper. By default set to 2g.", shortName="varMemory", required=false)
  var variantCallerMemory = 2

  @Argument(doc="expand each target in input intervals by the specified number of bases. By default set to 50 bases.", shortName="expand", required=false)
  var expandIntervals = 50

  @Argument(doc="pipeline memory limit. By default set to 2g.", shortName="pipeMemory", required=false)
  var pipelineMemoryLimit = 2

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.intervals = List(qscript.pipeline.getProject.getIntervalList)
    this.memoryLimit = pipelineMemoryLimit
  }

  def script() {
    pipeline = PicardPipeline.parse(qscript.yamlFile)

    val projectBase = qscript.pipeline.getProject.getName
    val bamType = "cleaned"

    val writeBamList = new ListWriterFunction
    writeBamList.inputFiles = qscript.pipeline.getSamples.filter(_.getBamFiles.contains(bamType)).map(_.getBamFiles.get(bamType)).toList
    writeBamList.listFile = projectBase +".bam.list"
    writeBamList.jobOutputFile = writeBamList.listFile + ".out"
    add(writeBamList)

    val flankIntervals = projectBase + ".flanks.intervals"

    if (qscript.expandIntervals > 0) {
      val ei = new ExpandIntervals(
        qscript.pipeline.getProject.getIntervalList,
        1,
        qscript.expandIntervals,
        flankIntervals,
        qscript.pipeline.getProject.getReferenceFile,
        "INTERVALS",
        "INTERVALS")
      ei.jobOutputFile = ei.outList + ".out"

      add(ei)
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0)
        this.intervals :+= flankIntervals
    }

    val call = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    call.input_file = List(writeBamList.listFile)
    call.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
    call.downsample_to_coverage = 600
    call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    call.GSA_PRODUCTION_ONLY = true
    call.out = projectBase + ".unfiltered.vcf"
    call.jobOutputFile = call.out + ".out"
    call.scatterCount = qscript.variantCallerScatterCount
    call.memoryLimit = qscript.variantCallerMemory
    add(call)

    val selectSNPs = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectSNPs.selectSNPs = true
    selectSNPs.rodBind :+= RodBind("variant", "VCF", call.out)
    selectSNPs.out = projectBase + ".snps.unfiltered.vcf"
    selectSNPs.jobOutputFile = selectSNPs.out + ".out"
    add(selectSNPs)

    val selectIndels = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectIndels.selectIndels = true
    selectIndels.rodBind :+= RodBind("variant", "VCF", call.out)
    selectIndels.out = projectBase + ".indels.unfiltered.vcf"
    selectIndels.jobOutputFile = selectIndels.out + ".out"
    add(selectIndels)

    val filterSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterSNPs.variantVCF = selectSNPs.out
    filterSNPs.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    filterSNPs.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    filterSNPs.clusterWindowSize = 10
    filterSNPs.clusterSize = 3
    filterSNPs.out = projectBase + ".snps.filtered.vcf"
    filterSNPs.jobOutputFile = filterSNPs.out + ".out"
    add(filterSNPs)

    val filterIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterIndels.variantVCF = selectIndels.out
    filterIndels.filterName = List("Indel_QUAL", "Indel_SB", "Indel_QD")
    filterIndels.filterExpression = List("\"QUAL<30.0\"", "\"SB>-1.0\"", "\"QD<2.0\"")
    filterIndels.out = projectBase + ".indels.filtered.vcf"
    filterIndels.jobOutputFile = filterIndels.out + ".out"
    add(filterIndels)

    val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineSNPsIndels.rodBind :+= RodBind("indels", "VCF", filterIndels.out)
    combineSNPsIndels.rodBind :+= RodBind("snps", "VCF", filterSNPs.out)
    combineSNPsIndels.rod_priority_list = "indels,snps"
    combineSNPsIndels.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineSNPsIndels.assumeIdenticalSamples = true
    combineSNPsIndels.out = projectBase + ".unannotated.vcf"
    combineSNPsIndels.jobOutputFile = combineSNPsIndels.out + ".out"
    add(combineSNPsIndels)

    val annotate = new GenomicAnnotator with CommandLineGATKArgs with ExpandedIntervals
    annotate.rodBind :+= RodBind("variant", "VCF", combineSNPsIndels.out)
    annotate.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.pipeline.getProject.getRefseqTable)
    annotate.rodToIntervalTrackName = "variant"
    annotate.out = projectBase + ".vcf"
    annotate.jobOutputFile = annotate.out + ".out"
    add(annotate)

    val targetEval = new VariantEval with CommandLineGATKArgs
    targetEval.rodBind :+= RodBind("eval", "VCF", annotate.out)
    targetEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
    targetEval.doNotUseAllStandardStratifications = true
    targetEval.doNotUseAllStandardModules = true
    targetEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    targetEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
    targetEval.out = projectBase + ".eval"
    targetEval.jobOutputFile = targetEval.out + ".out"
    add(targetEval)

    if (qscript.expandIntervals > 0) {
      val flanksEval = new VariantEval with CommandLineGATKArgs
      flanksEval.rodBind :+= RodBind("eval", "VCF", annotate.out)
      flanksEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
      flanksEval.intervals = List(flankIntervals)
      flanksEval.doNotUseAllStandardStratifications = true
      flanksEval.doNotUseAllStandardModules = true
      flanksEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
      flanksEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
      flanksEval.out = projectBase + ".flanks.eval"
      flanksEval.jobOutputFile = flanksEval.out + ".out"
      add(flanksEval)
    }
  }
}
