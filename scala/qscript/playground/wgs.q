import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class WGSpipeline extends QScript {
  qscript =>

  @Input(doc="File with list of Bams", shortName = "bams", required = true)
  var bamList: File = _

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="Project Name", shortName = "P", required = true)
  var project: String =_

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=true)
  var outputTmpDir: String = _

  @Input(doc="reference file", shortName = "R", required = true)
  var reference: File = _

  @Input(doc="dbSNP", shortName = "D", required = true)
  var dbSNP: File = _

  @Input(doc="gsa-pipeline directory key", shortName =  "dirKey", required = true)
  var dirkey: String =_

  @Input(doc="version id in format ### (ie 001)", shortName = "v", required = true)
  var version: String =_

  @Input(doc="Flag for running the entire genome. Otherwise a QC run of chr 20 is all that is run.", shortName = "freeze", required = false)
  var freeze: Boolean = false


  private val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf"

  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
    val g1k = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf"
    val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf"

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(3)
  }

  def script = {
    var chr: Int = 20
    var runType: String = "QC"
    if (qscript.freeze) {
      var chr: Int =20
      var runType = "freeze"
    }  //TODO: Impliment whole genome once infrastructure for waiting for tmp space is available.
    // for now this will only run on chromosome 20 because otherwise we might run  into space errors.
    else{
      var chr: Int = 20
      var runType = "QC"
    }
    val dirbase = qscript.dirkey + "/" + runType + "_v" +qscript.version
    val base = dirbase + "/Chunks/" + qscript.project + "." + runType
    val projectBase = dirbase + "/" + qscript.project + "." + runType
    val basesPerJob: Int = 3000000
    val lastBase: Int = 63025520
    var start: Int = 1
    var stop: Int = start - 1 + basesPerJob
    if( stop > lastBase ) { stop = lastBase }
    var jobNumber: Int = 1
    var mergeList: List[List[_]] = Nil
    while( jobNumber < (lastBase.toFloat / basesPerJob.toFloat) + 1.0) {
      var interval = "%d:%d-%d".format(chr, start, stop)
      val baseTmpName: String = qscript.outputTmpDir + "/" +base  + "chunk_" +jobNumber
      val baseJobName: String = qscript.project + "_chunk_" + jobNumber
      // 1.) Create cleaning targets
      val target = new RealignerTargetCreator with CommandLineGATKArgs
      target.memoryLimit = Some(3)
      target.input_file :+= qscript.bamList
      target.intervalsString :+= interval
      target.out = new File(baseTmpName + ".target.intervals")
      target.mismatchFraction = Some(0.0)
      target.maxIntervalSize = Some(700)
      target.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      target.jobName = baseJobName + ".target"
      target.jobOutputFile = new File(baseTmpName + ".target.out")
      target.isIntermediate = true
      // 2.) Clean without SW
      val clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + ".cleaned.bam")
      clean.memoryLimit = Some(4)
      clean.input_file :+= qscript.bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = target.out
      clean.jobOutputFile = new File (baseTmpName + ".clean.out"  )
      clean.out = cleanedBam
      clean.doNotUseSW = true
      clean.simplifyBAM = true
      clean.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      clean.jobName = baseJobName + ".clean"
      clean.isIntermediate = true
      clean.compress = Some(0)
      clean.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )

      val call = new UnifiedGenotyper with CommandLineGATKArgs
      call.out = new File("/humgen/gsa-pipeline/"+base+ "chunk_" +jobNumber + ".vcf")
      call.dcov = Some( 50 )
      call.stand_call_conf = Some( 4.0 )
      call.stand_emit_conf = Some( 4.0 )
      call.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
      call.jobName =  baseJobName  +".call"
      call.exactCalculation = org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel.ExactCalculation.LINEAR_EXPERIMENTAL
      call.jobOutputFile = new File ("/humgen/gsa-pipeline/"+base+ "chunk_" +jobNumber  + ".call.out"    )
      call.input_file = List(clean.out)
      call.intervalsString :+= interval
      call.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )

      mergeList:+= List(baseJobName, call.out)
      add(target, clean, call)

      start += basesPerJob
      stop += basesPerJob
      if( stop > lastBase ) { stop = lastBase }
      jobNumber += 1
    }

    val combineVCFs = new CombineVariants with CommandLineGATKArgs
    combineVCFs.priority = ""
    for (c <- mergeList ){
    combineVCFs.rodBind :+= RodBind(c(0).toString, "VCF", c(1).toString)
    combineVCFs.priority += c(0).toString +","
    }
    combineVCFs.variantMergeOptions = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION
    combineVCFs.jobName = qscript.project + ".combine"
    combineVCFs.out = new File("/humgen/gsa-pipeline/" + projectBase + ".merged.vcf")
    combineVCFs.jobOutputFile = new File ("/humgen/gsa-pipeline/" +projectBase + ".combine.out")
    combineVCFs.genotypemergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.PRIORITIZE
    combineVCFs.intervalsString :+= "20:1-63025520"   ///TODO impliment only using the chr 20 or nothing if freeze



    val sv = new SelectVariants with CommandLineGATKArgs
    sv.selectIndels = true
    sv.out = swapExt(combineVCFs.out, "merged.vcf", "indelsOnly.vcf")
    sv.rodBind :+= RodBind("variant", "VCF", combineVCFs.out)
    sv.jobName = qscript.project + ".SV"
    sv.jobOutputFile = swapExt(combineVCFs.jobOutputFile, "combine.out", "sv.out")

    val filter = new VariantFiltration with CommandLineGATKArgs
      filter.intervalsString :+= "20:1-63025520"                 ///TODO impliment only using the chr 20 or nothing if freeze
      filter.variantVCF = sv.out
      filter.out = swapExt(sv.out, ".vcf", ".filtered.vcf")
      filter.filterName ++= List("HARD_TO_VALIDATE")            //ToDO check to make sure that this is what Guillermo recomends
      filter.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
      filter.jobName = qscript.project + ".VF"
    filter.jobOutputFile = swapExt(sv.jobOutputFile, ".sv.out", ".filter.out")


    val recombine = new CombineVariants with CommandLineGATKArgs
    recombine.rodBind :+= RodBind("indels", "VCF", sv.out)
    recombine.rodBind :+= RodBind("all", "VCF", combineVCFs.out)
    recombine.priority = "indels, all"
    recombine.genotypemergeoption = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.PRIORITIZE
    recombine.out = swapExt(combineVCFs.out, "vcf", "filtered.vcf")
    recombine.jobName = qscript.project + ".recombine"
    recombine.jobOutputFile = swapExt(combineVCFs.jobOutputFile, "combine.out", "recombine.out")

    val cr = new ContrastiveRecalibrator with CommandLineGATKArgs
      cr.rodBind :+= RodBind("hapmap", "VCF", qscript.hapmap)
      cr.rodBind :+= RodBind("1kg", "VCF", qscript.omni)
      cr.rodBind :+= RodBind("input", "VCF", recombine.out)
      cr.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
      cr.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      cr.jobName = qscript.project + ".cr"
      cr.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".cr.out")
      cr.intervalsString :+= "20"                 ///TODO impliment only using the chr 20 or nothing if freeze
      cr.allPoly = true
      cr.tranches_file  =swapExt(filter.out, "merged.filtered.vcf", "tranches")
      cr.recal_file = swapExt(filter.out, "merged.filtered.vcf", "")
      cr.tranche ++= List("0.1","0.5","0.75","1.0","1.25","1.35","1.5","1.7","1.8","1.9","2.0","2.1","2.2","2.3","2.5","3.0","5.0","8.0","10.0")

    val ar = new ApplyRecalibration with CommandLineGATKArgs
    ar.tranches_file = cr.tranches_file
    ar.recal_file = cr.recal_file
    ar.ignore_filter++= List("HARD_TO_VALIDATE")
    ar.out = swapExt(filter.out, "vcf", "recalibrated.vcf")
      ar.jobName = qscript.project + ".ar"
      ar.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".ar.out")

    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.jobName = qscript.project+".eval"
    stdEval.jobOutputFile = swapExt(combineVCFs.jobOutputFile, ".combine.out", ".eval.out")
    stdEval.noST = true
    stdEval.noEV = true
    stdEval.evalModule ++= List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    stdEval.stratificationModule ++= List("EvalRod", "CompRod", "Novelty")
    stdEval.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
    stdEval.rodBind :+= RodBind("eval", "VCF", ar.out)
    stdEval.out = swapExt(ar.out, ".vcf", ".eval")

    add(combineVCFs, filter, sv, recombine, cr, ar, stdEval)
  }
}
