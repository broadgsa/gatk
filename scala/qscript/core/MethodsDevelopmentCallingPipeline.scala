import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport


  // ToDos:
  // reduce the scope of the datasets so the script is more nimble
  // create gold standard BAQ'd bam files, no reason to always do it on the fly

  // Analysis to add at the end of the script:
  // auto generation of the cluster plots
  // spike in NA12878 to the exomes and to the lowpass, analysis of how much of her variants are being recovered compared to single sample exome or HiSeq calls
  // produce Kiran's Venn plots based on comparison between new VCF and gold standard produced VCF


class MethodsDevelopmentCallingPipeline extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=true)
  var outputDir: String = "./"

  @Argument(shortName="skipCalling", doc="skip the calling part of the pipeline and only run VQSR on preset, gold standard VCF files", required=false)
  var skipCalling: Boolean = false

  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil

  @Argument(shortName="skipGoldStandard", doc="doesn't run the pipeline with the goldstandard VCF files for comparison", required=false)
  var skipGoldStandard: Boolean = false

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="eval", doc="adds the VariantEval walker to the pipeline", required=false)
  var eval: Boolean = false

  @Argument(shortName="indels", doc="calls indels with the Unified Genotyper", required=false)
  var callIndels: Boolean = false

  @Argument(shortName="LOCAL_ET", doc="Doesn't use the AWS S3 storage for ET option", required=false)
  var LOCAL_ET: Boolean = false

  @Argument(shortName="mbq", doc="The minimum Phred-Scaled quality score threshold to be considered a good base.", required=false)
  var minimumBaseQuality: Int = -1

  @Argument(shortName="deletions", doc="Maximum deletion fraction allowed at a site to call a genotype.", required=false)
  var deletions: Double = -1

  @Argument(shortName="sample", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = Nil



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
          val trancheTarget: Double,
          val isLowpass: Boolean) {
    val name = qscript.outputDir + baseName
    val clusterFile = new File(name + ".clusters")
    val rawVCF = new File(name + ".raw.vcf")
    val rawIndelVCF = new File(name + ".raw.indel.vcf")
    val filteredIndelVCF = new File(name + ".filtered.indel.vcf")
    val recalibratedVCF = new File(name + ".recalibrated.vcf")
    val tranchesFile = new File(name + ".tranches")
    val vqsrRscript = name + ".vqsr.r"
    val recalFile = new File(name + ".tranches.recal")
    val goldStandardRecalibratedVCF = new File(name + "goldStandard.recalibrated.vcf")
    val goldStandardTranchesFile = new File(name + "goldStandard.tranches")
    val goldStandardRecalFile = new File(name + "goldStandard.tranches.recal")
    val evalFile = new File(name + ".snp.eval")
    val evalIndelFile = new File(name + ".indel.eval")
    val goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    val goldStandardClusterFile = new File(goldStandardName + ".clusters")
  }

  val hg19 = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")  
  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val dbSNP_hg18_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_hg18.rod"            // Special case for NA12878 collections that can't use 132 because they are part of it.
  val dbSNP_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b36.rod"
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"
  val dbSNP_b37_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf"              // Special case for NA12878 collections that can't use 132 because they are part of it.
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val training_hapmap_b37 = "/humgen/1kg/processing/pipeline_test_bams/hapmap3.3_training_chr20.vcf"
  val omni_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b36.vcf"
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val indelMask_b36 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b36.bed"
  val indelMask_b37 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b37.bed"

  val lowPass: Boolean = true
  val indels: Boolean = true

  val queueLogDir = ".qlog/"

  val targetDataSets: Map[String, Target] = Map(
    "HiSeq" -> new Target("NA12878.HiSeq", hg18, dbSNP_hg18_129, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals", 2.07, 99.0, !lowPass),
    "HiSeq19" -> new Target("NA12878.HiSeq19", hg19, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/hiseq19/analysis/snps/NA12878.HiSeq19.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.3, 99.0, !lowPass),
    "GA2hg19" -> new Target("NA12878.GA2.hg19", hg19, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.GA2.WGS.bwa.cleaned.hg19.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/hiseq19/analysis/snps/NA12878.GA2.hg19.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.3, 99.0, !lowPass),
    "WEx" -> new Target("NA12878.WEx", hg18, dbSNP_hg18_129, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, 97.0, !lowPass),
    "WExTrio" -> new Target("CEUTrio.WEx", hg19, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.bwa.cleaned.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, 97.0, !lowPass),
    "WGSTrio" -> new Target("CEUTrio.WGS", hg19, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.3, 99.0, !lowPass),
    "FIN" -> new Target("FIN", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, 99.0, lowPass),
    "TGPWExGdA" -> new Target("1000G.WEx.GdA", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"),        // BUGBUG: reduce from 60 to 20 people
              new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, 99.0, !lowPass),
    "LowPassN60" -> new Target("lowpass.N60", b36, dbSNP_b36, hapmap_b36, indelMask_b36,
              new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
              new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, 99.0, lowPass),          // chunked interval list to use with Queue's scatter/gather functionality
    "LowPassEUR363Nov" -> new Target("EUR.nov2010", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, 99.0, lowPass)
  )


  def script = {

    // Selects the datasets in the -dataset argument and adds them to targets.
    var targets: List[Target] = List()
    if (!datasets.isEmpty)
      for (ds <- datasets)
        targets ::= targetDataSets(ds)
    else                                                // If -dataset is not specified, all datasets are used.
      for (targetDS <- targetDataSets.valuesIterator)
        targets ::= targetDS

    val goldStandard = true
    for (target <- targets) {
      if( !skipCalling ) {
        if (callIndels) add(new indelCall(target), new indelFilter(target), new indelEvaluation(target))
        add(new snpCall(target))
        add(new VQSR(target, !goldStandard))
        add(new applyVQSR(target, !goldStandard))
        if (eval) add(new snpEvaluation(target))
      }
      if ( !skipGoldStandard ) {
        add(new VQSR(target, goldStandard))
        add(new applyVQSR(target, goldStandard))
      }
    }
  }


  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    jarFile = gatkJarFile;
    memoryLimit = 4;
    phone_home = if ( LOCAL_ET ) GATKRunReport.PhoneHomeOption.STANDARD else GATKRunReport.PhoneHomeOption.AWS_S3
  }

  def bai(bam: File) = new File(bam + ".bai")
  val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

  // 1.) Unified Genotyper Base
  class GenotyperBase (t: Target) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
    this.dcov = if ( t.isLowpass ) { 50 } else { 250 }
    this.stand_call_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.stand_emit_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.input_file :+= t.bamList
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }

  // 1a.) Call SNPs with UG
  class snpCall (t: Target) extends GenotyperBase(t) {
    if (minimumBaseQuality >= 0)
      this.min_base_quality_score = minimumBaseQuality
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = t.rawVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY}
    this.analysisName = t.name + "_UGs"
    this.jobName =  queueLogDir + t.name + ".snpcall"
    this.A ++= List("FisherStrand")
  }

  // 1b.) Call Indels with UG
  class indelCall (t: Target) extends GenotyperBase(t) {
    this.out = t.rawIndelVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    this.analysisName = t.name + "_UGi"
    this.jobName =  queueLogDir + t.name + ".indelcall"
  }

  // 2.) Hard Filtering for indels
  class indelFilter (t: Target) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 2
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 10
    this.filterName ++= List("HARD_TO_VALIDATE")
    this.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    this.variantVCF = t.rawIndelVCF
    this.out = t.filteredIndelVCF
    this.filterName ++= List("LowQual", "StrandBias", "QualByDepth", "HomopolymerRun")
    if (t.isLowpass)
      this.filterExpression ++= List("\"QUAL<30.0\"", "\"SB>=-1.0\"", "\"QD<1.0\"", "\"HRun>=15\"")
    else
      this.filterExpression ++= List("\"QUAL<50.0\"", "\"SB>=-1.0\"", "\"QD<5.0\"", "\"HRun>=15\"")
    this.analysisName = t.name + "_VF"
    this.jobName =  queueLogDir + t.name + ".indelfilter"
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(t: Target, goldStandard: Boolean) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.rawVCF } )
    this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile, "known=false,training=true,truth=true,prior=15.0")
    if( t.hapmapFile.contains("b37") )
      this.rodBind :+= RodBind("omni", "VCF", omni_b37, "known=false,training=true,truth=true,prior=12.0")
    else if( t.hapmapFile.contains("b36") )
      this.rodBind :+= RodBind("omni", "VCF", omni_b36, "known=false,training=true,truth=true,prior=12.0")
    if (t.dbsnpFile.endsWith(".rod"))
      this.rodBind :+= RodBind("dbsnp", "DBSNP", t.dbsnpFile, "known=true,training=false,truth=false,prior=10.0")
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile, "known=true,training=false,truth=false,prior=10.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "HRun", "FS")
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile }
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.rscript_file = t.vqsrRscript
    this.analysisName = t.name + "_VQSR"
    this.jobName =  queueLogDir + t.name + ".VQSR"
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  class applyVQSR (t: Target, goldStandard: Boolean) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.rawVCF } )
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile}
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.ts_filter_level = t.trancheTarget
    this.out = t.recalibratedVCF
    this.analysisName = t.name + "_AVQSR"
    this.jobName =  queueLogDir + t.name + ".applyVQSR"
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class EvalBase(t: Target) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.reference_sequence = t.reference
    this.rodBind :+= RodBind("comphapmap", "VCF", t.hapmapFile)
    this.intervalsString ++= List(t.intervals)
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
    this.sample = samples
  }

  // 5a.) SNP Evaluation (OPTIONAL) based on the cut vcf
  class snpEvaluation(t: Target) extends EvalBase(t) {
    if (t.reference == b37 || t.reference == hg19) this.rodBind :+= RodBind("compomni", "VCF", omni_b37)
    this.rodBind :+= RodBind("eval", "VCF", t.recalibratedVCF )
    this.out =  t.evalFile
    this.analysisName = t.name + "_VEs"
    this.jobName =  queueLogDir + t.name + ".snp.eval"
  }

  // 5b.) Indel Evaluation (OPTIONAL)
  class indelEvaluation(t: Target) extends EvalBase(t) {
    this.rodBind :+= RodBind("eval", "VCF", t.filteredIndelVCF)
    this.evalModule :+= "IndelStatistics"
    this.out =  t.evalIndelFile
    this.analysisName = t.name + "_VEi"
    this.jobName =  queueLogDir + queueLogDir + t.name + ".indel.eval"
  }
}
