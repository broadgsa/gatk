import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;

class MethodsDevelopmentCallingPipeline extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=true)
  var outputDir: String = "./"

  @Argument(shortName="skipCalling", doc="If true, skip the calling part of the pipeline and only run VQSR on preset, gold standard VCF files", required=false)
  var skipCalling: Boolean = false

  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil

  @Argument(shortName="skipGoldStandard", doc="runs the pipeline with the goldstandard VCF files for comparison", required=false)
  var skipGoldStandard: Boolean = false

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="noMASK", doc="turns off MASK calculation", required=false)
  var noMASK: Boolean = false

  @Argument(shortName="eval", doc="adds the VariantEval walker to the pipeline", required=false)
  var eval: Boolean = false


  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "INFO"; jarFile = gatkJarFile; memoryLimit = Some(3); }

  class Target(
          val baseName: String,
          val reference: File,
          val dbsnpFile: String,
          val hapmapFile: String,
          val maskFile: String,
          val bamList: File,
          val goldStandard_VCF: File,
          val intervals: String,
          val titvTarget: Double,
          val isLowpass: Boolean) {
    val name = qscript.outputDir + baseName
    val clusterFile = new File(name + ".clusters")
    val rawVCF = new File(name + ".raw.vcf")
    val filteredVCF = new File(name + ".filtered.vcf")
    val titvRecalibratedVCF = new File(name + ".titv.recalibrated.vcf")
    val tsRecalibratedVCF = new File(name + ".ts.recalibrated.vcf")
    val goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    val goldStandardClusterFile = new File(goldStandardName + ".clusters")
  }

  val hg19 = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")  
  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val dbSNP_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_130_hg18.rod"
  val dbSNP_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_130_b36.rod"
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val indelMask_b36 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b36.bed"
  val indelMask_b37 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b37.bed"

  // ToDos:
  // reduce the scope of the datasets so the script is more nimble
  // figure out how to give names to all the Queue-LSF logs (other than Q-1931@node1434-24.out) so that it is easier to find logs for certain steps
  // create gold standard BAQ'd bam files, no reason to always do it on the fly

  // Analysis to add at the end of the script:
  // auto generation of the cluster plots
  // spike in NA12878 to the exomes and to the lowpass, analysis of how much of her variants are being recovered compared to single sample exome or HiSeq calls
  // produce Kiran's Venn plots based on comparison between new VCF and gold standard produced VCF

  val lowPass: Boolean = true

  val targetDataSets: Map[String, Target] = Map(
    "HiSeq" -> new Target("NA12878.HiSeq", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals", 2.07, !lowPass),
    "FIN" -> new Target("FIN", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass),
    "WEx" -> new Target("NA12878.WEx", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, !lowPass),
    "TGPWExGdA" -> new Target("1000G.WEx.GdA", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"),        // BUGBUG: reduce from 60 to 20 people
              new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass),
    "LowPassN60" -> new Target("lowpass.N60", b36, dbSNP_b36, hapmap_b36, indelMask_b36,
              new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
              new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, lowPass),          // chunked interval list to use with Queue's scatter/gather functionality
    "LowPassAugust" -> new Target("ALL.august.v4", b37, dbSNP_b37, hapmap_b37, indelMask_b37,                               // BUGBUG: kill this, it is too large
              new File("/humgen/1kg/processing/allPopulations_chr20_august_release.cleaned.merged.bams/ALL.cleaned.merged.list"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass),
    "LowPassEUR363Nov" -> new Target("EUR.nov2010", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass),
    "WExTrio" -> new Target("NA12878Trio.WEx", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
        new File("/humgen/gsa-scr1/carneiro/prj/trio/data/NA12878Trio.WEx.hg19.recal.bam"),
        new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass)
  )


  def script = {

    // Selects the datasets in the -dataset argument and adds them to targets.
    var targets: List[Target] = List()
    if (!datasets.isEmpty)
      for (ds <- datasets)
        targets ::= targetDataSets(ds)                  // Could check if ds was mispelled, but this way an exception will be thrown, maybe it's better this way?
    else                                                // If -dataset is not specified, all datasets are used.
      for (targetDS <- targetDataSets.valuesIterator)   // for Scala 2.7 or older, use targetDataSets.values
        targets ::= targetDS

    val goldStandard = true
    for (target <- targets) {
      if( !skipCalling ) {
        add(new UnifiedGenotyper(target))
        add(new VariantFiltration(target))
        add(new GenerateVariantClusters(target, !goldStandard))
        add(new VariantRecalibratorTiTv(target, !goldStandard))
        add(new VariantRecalibratorNRS(target, !goldStandard))
        if (eval) add(new VariantEvaluation(target))
      }
      if ( !skipGoldStandard ) {
        add(new GenerateVariantClusters(target, goldStandard))
        add(new VariantRecalibratorTiTv(target, goldStandard))
        add(new VariantRecalibratorNRS(target, goldStandard))
      }
    }
  }

  def bai(bam: File) = new File(bam + ".bai")

  val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

  // 1.) Call SNPs with UG
  class UnifiedGenotyper(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
    this.dcov = Some( if ( t.isLowpass ) { 50 } else { 250 } )
    this.stand_call_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.stand_emit_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.input_file :+= t.bamList
    this.out = t.rawVCF
    this.baq = Some( if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE})
    this.analysisName = t.name + "_UG"
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }

  // 2.) Filter SNPs
  class VariantFiltration(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 10
    this.variantVCF = t.rawVCF
    this.out = t.filteredVCF
    this.filterName ++= List("HARD_TO_VALIDATE")
    this.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    this.analysisName = t.name + "_VF"
    if (!noMASK) {
      this.rodBind :+= RodBind("mask", "Bed", t.maskFile)
      this.maskName = "InDel"
    }
  }

  // 3.) VQSR part1 Generate Gaussian clusters based on truth sites
  class GenerateVariantClusters(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.GenerateVariantClusters with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      if( t.hapmapFile.contains("b37") )
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf")
      this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      this.analysisName = name + "_GVC"
      this.intervalsString ++= List(t.intervals)
      this.qual = Some(350) // clustering parameters to be updated soon pending new experimentation results
      this.std = Some(3.5)
      this.mG = Some(10)
      this.ignoreFilter ++= FiltersToIgnore
      if (t.dbsnpFile.endsWith(".rod"))
        this.DBSNP = new File(t.dbsnpFile)
      else if (t.dbsnpFile.endsWith(".vcf"))
        this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }

  // 4.) VQSR part2 Calculate new LOD for all input SNPs by evaluating the Gaussian clusters
  class VariantRecalibratorBase(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.VariantRecalibrator with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      if( t.hapmapFile.contains("b37") )
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf")
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("truth", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      this.analysisName = name + "_VR"
      this.intervalsString ++= List(t.intervals)
      this.ignoreFilter ++= FiltersToIgnore
      this.ignoreFilter ++= List("HARD_TO_VALIDATE")
      this.target_titv = Some(t.titvTarget)
      if (t.dbsnpFile.endsWith(".rod"))
        this.DBSNP = new File(t.dbsnpFile)
      else if (t.dbsnpFile.endsWith(".vcf"))
        this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }

  // 4a.) Choose VQSR tranches based on novel ti/tv
  class VariantRecalibratorTiTv(t: Target, goldStandard: Boolean) extends VariantRecalibratorBase(t, goldStandard) {
      this.tranche ++= List("0.1", "1.0", "10.0", "100.0")
      this.out = new File(this.name + ".titv.recalibrated.vcf")
      this.tranchesFile = new File(this.name + ".titv.tranches")
  }

  // 4b.) Choose VQSR tranches based on sensitivity to truth set
  class VariantRecalibratorNRS(t: Target, goldStandard: Boolean) extends VariantRecalibratorBase(t, goldStandard) {
      this.sm = Some(org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY)
      this.tranche ++= List("0.1", "1.0", "10.0", "100.0")
      this.out = new File(this.name + ".ts.recalibrated.vcf")
      this.priorDBSNP = Some(2.0)
      this.priorHapMap = Some(2.0)
      this.prior1KG = Some(2.0)    
      this.tranchesFile = new File(this.name + ".ts.tranches")
  }

  // 5.) Variant Evaluation (OPTIONAL!) based on the sensitivity recalibrated vcf
  class VariantEvaluation(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantEval with UNIVERSAL_GATK_ARGS {
      val name: String = t.name
      this.reference_sequence = t.reference
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("eval", "VCF", t.tsRecalibratedVCF)
      this.analysisName = name + "_VR"
      this.intervalsString ++= List(t.intervals)
      this.reportType = Some(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.R)
      this.reportLocation = new File(t.name + ".eval")
      this.noStandard = true
      this.evalModule ++= List("TiTvVariantEvaluator", "CountVariants", "GenotypeConcordance")
      if (t.dbsnpFile.endsWith(".rod"))
        this.DBSNP = new File(t.dbsnpFile)
	    else if (t.dbsnpFile.endsWith(".vcf"))
        this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }
}
