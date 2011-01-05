import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
//import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
//import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
//import org.apache.commons.io.FilenameUtils;

class tdPipeline extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=true)
  var outputDir: String = "./"

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var useBAQ: Boolean = true

  @Argument(shortName="noMask", doc="turns off MASK calculation", required=false)
  var useMask: Boolean = true


  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "INFO"; jarFile = gatkJarFile; memoryLimit = Some(4);}

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
    val goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    val goldStandardClusterFile = new File(goldStandardName + ".clusters")
  }

  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  // Define the target datasets here
  val  lowPass: Boolean = true
  val HiSeq = new Target("NA12878.HiSeq",
        hg18,
        "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_hg18.rod",
        "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.hg18_fwd.vcf",
        "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.bed",
        new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
        new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals",
        2.07,
        !lowPass)

  val LowPassFIN79Nov = new Target("FIN.nov2010",
        b37,
        "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf",
        "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf",
        "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.hg18.bed",
        new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
        new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals",
        2.3,
        lowPass)

  /*
   * These sources need to be updated, never used.

  val TGPWExGdA = new Target("1000G.WEx.GdA", b37, "b37",
        new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"), // BUGBUG: reduce from 60 to 20 people
        new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass)
  val WEx = new Target("NA12878.WEx", hg18, "hg18",
        new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
        new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, !lowPass)
  val LowPassN60 = new Target("lowpass.N60", b36, "b36",                                                              // which reference the data is aligned to
        new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
        new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, lowPass)           // chunked interval list to use with Queue's scatter/gather functionality
  val LowPassAugust = new Target("ALL.august.v4", b37, "b37", // BUGBUG: kill this, it is too large
        new File("/humgen/1kg/processing/allPopulations_chr20_august_release.cleaned.merged.bams/ALL.cleaned.merged.list"),
        new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass)
  val LowPassEUR363Nov = new Target("EUR.nov2010", b37, "b37",
        new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
        new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass)
  */
  val targets = List(HiSeq, LowPassFIN79Nov)

  def script = {
      val goldStandard = true
      for (target <- targets) {
        add(new UnifiedGenotyper(target))
        add(new VariantFiltration(target))
        add(new GenerateVariantClusters(target, !goldStandard))
        add(new VariantRecalibratorTiTv(target, !goldStandard))
        add(new VariantRecalibratorNRS(target, !goldStandard))
        add(new VariantEvaluation(target))
      }
  }

  val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

  // 1.) Call SNPs with UG
  class UnifiedGenotyper(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    if (t.dbsnpFile.endsWith(".rod")) this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf")) this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
    this.jobQueue = "hour"
    this.dcov = Some( if ( t.isLowpass ) { 50 } else { 250 } )
    this.input_file :+= t.bamList
    this.out = t.rawVCF
    if  (useBAQ)
      this.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE)
    else
      this.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF)
    this.stand_call_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.stand_emit_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.analysisName = t.name + "_UG"
  }

  // 2.) Filter SNPs
  class VariantFiltration(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 10
    this.jobQueue = "hour"
    this.variantVCF = t.rawVCF
    this.out = t.filteredVCF
    if (useMask) {
      this.rodBind :+= RodBind("mask", "Bed", t.maskFile)
      this.maskName = "InDel"
    }
    this.filterName ++= List("HARD_TO_VALIDATE")
    this.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    this.analysisName = t.name + "_VF"
  }

  // 3.) VQSR part1 Generate Gaussian clusters based on truth sites
  class GenerateVariantClusters(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.GenerateVariantClusters with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      if (t.dbsnpFile.endsWith(".rod")) this.DBSNP = new File(t.dbsnpFile)
	    else if (t.dbsnpFile.endsWith(".vcf")) this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      this.analysisName = name + "_GVC"
      this.intervalsString ++= List(t.intervals)
      this.qual = Some(300) // clustering parameters to be updated soon pending new experimentation results
      this.std = Some(3.5)
      this.mG = Some(16)
      this.ignoreFilter ++= FiltersToIgnore
  }

  // 4.) VQSR part2 Calculate new LOD for all input SNPs by evaluating the Gaussian clusters
  class VariantRecalibratorBase(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.VariantRecalibrator with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      if (t.dbsnpFile.endsWith(".rod")) this.DBSNP = new File(t.dbsnpFile)
	    else if (t.dbsnpFile.endsWith(".vcf")) this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("truth", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      this.analysisName = name + "_VR"
      this.intervalsString ++= List(t.intervals)
      this.ignoreFilter ++= FiltersToIgnore
      this.ignoreFilter ++= List("HARD_TO_VALIDATE")
      this.priorDBSNP = Some(2.0)
      this.priorHapMap = Some(2.0)
      this.target_titv = t.titvTarget
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
      this.tranchesFile = new File(this.name + ".ts.tranches")
  }

  // 5.) Variant Evaluation generates the csv files for Genome Concordance analysis
  class VariantEvaluation(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantEval with UNIVERSAL_GATK_ARGS {
      val name: String = t.name
      this.reference_sequence = t.reference
      if (t.dbsnpFile.endsWith(".rod")) this.DBSNP = new File(t.dbsnpFile)
	    else if (t.dbsnpFile.endsWith(".vcf")) this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
      this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("input", "VCF", t.filteredVCF)
      this.analysisName = name + "_VR"
      this.intervalsString ++= List(t.intervals)
      this.reportType = Some(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.R)
      this.reportLocation = new File(t.name + ".eval")
  }
}
