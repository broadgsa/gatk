package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport

class MethodsDevelopmentCallingPipeline extends QScript {
  qscript =>

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./"

  @Argument(shortName="skipCalling", doc="skip the calling part of the pipeline and only run VQSR on preset, gold standard VCF files", required=false)
  var skipCalling: Boolean = false

  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil

  @Argument(shortName="runGoldStandard", doc="run the pipeline with the goldstandard VCF files for comparison", required=false)
  var runGoldStandard: Boolean = false

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="noIndels", doc="do not call indels with the Unified Genotyper", required=false)
  var noIndels: Boolean = false

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
          val isLowpass: Boolean,
          val isExome: Boolean,
          val nSamples: Int) {
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

  val b37_decoy = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val hg19 = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val dbSNP_hg18_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_hg18.rod"            // Special case for NA12878 collections that can't use 132 because they are part of it.
  val dbSNP_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132.b36.excluding_sites_after_129.vcf"
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf"              // Special case for NA12878 collections that can't use 132 because they are part of it.
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val training_hapmap_b37 = "/humgen/1kg/processing/pipeline_test_bams/hapmap3.3_training_chr20.vcf"
  val omni_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b36.vcf"
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val indelMask_b36 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b36.bed"
  val indelMask_b37 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b37.bed"
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val badSites_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"

  val lowPass: Boolean = true
  val exome: Boolean = true
  val indels: Boolean = true

  val queueLogDir = ".qlog/"

  // BUGBUG: We no longer support b36/hg18 because several of the necessary files aren't available aligned to those references

  val targetDataSets: Map[String, Target] = Map(
    "NA12878_gold" -> new Target("NA12878.goldStandard", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/dev/carneiro/NA12878/data/goldStandard.list"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/NA12878/analysis/snps/NA12878.HiSeq19.filtered.vcf"),          // ** There is no gold standard for the gold standard **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.noChrY.hg19.intervals", 2.14, 99.0, lowPass, !exome, 391),
    "NA12878_wgs_b37" -> new Target("NA12878.HiSeq.WGS.b37", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/NA12878/analysis/snps/NA12878.HiSeq19.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.noChrY.hg19.intervals", 2.14, 99.0, !lowPass, !exome, 1),
    "NA12878_wgs_decoy" -> new Target("NA12878.HiSeq.WGS.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/NA12878/analysis/snps/NA12878.HiSeq19.filtered.vcf"),          // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.noChrY.hg19.intervals", 2.14, 99.0, !lowPass, !exome, 1),
    "NA12878_wgs_hg18" -> new Target("NA12878.HiSeq.WGS.hg18", hg18, dbSNP_hg18_129, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals", 2.14, 99.0, !lowPass, !exome, 1),
    "NA12878_wex_b37" -> new Target("NA12878.HiSeq.WEx.b37", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/seq/picard_aggregation/C339/NA12878/v3/NA12878.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 3.3, 98.0, !lowPass, exome, 1),
    "NA12878_wex_hg18" -> new Target("NA12878.HiSeq.WEx.hg18", hg18, dbSNP_hg18_129, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 3.3, 98.0, !lowPass, exome, 1),
    "NA12878_wex_decoy" -> new Target("NA12878.HiSeq.WEx.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 3.3, 98.0, !lowPass, exome, 1),
    "CEUTrio_wex_b37" -> new Target("CEUTrio.HiSeq.WEx.b37", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.bwa.cleaned.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 3.3, 98.0, !lowPass, exome, 3),
    "CEUTrio_wgs_b37" -> new Target("CEUTrio.HiSeq.WGS.b37", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.3, 99.0, !lowPass, !exome, 3),
    "CEUTrio_wex_decoy" -> new Target("CEUTrio.HiSeq.WEx.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.b37_decoy.list"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 3.3, 98.0, !lowPass, exome, 3),
    "CEUTrio_wgs_decoy" -> new Target("CEUTrio.HiSeq.WGS.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.list"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.3, 99.0, !lowPass, !exome, 3),
    "GA2hg19" -> new Target("NA12878.GA2.hg19", hg19, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.GA2.WGS.bwa.cleaned.hg19.bam"),
              new File("/humgen/gsa-hpprojects/dev/carneiro/NA12878/analysis/snps/NA12878.GA2.hg19.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 2.14, 99.0, !lowPass, !exome, 1),
    "FIN" -> new Target("FIN", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, 99.0, lowPass, !exome, 79),
    "TGPWExGdA" -> new Target("1000G.WEx.GdA", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"),        // BUGBUG: reduce from 60 to 20 people
              new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, 99.0, !lowPass, exome, 96),
    "LowPassN60" -> new Target("lowpass.N60", b36, dbSNP_b36, hapmap_b36, indelMask_b36,
              new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
              new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, 99.0, lowPass, !exome, 60),          // chunked interval list to use with Queue's scatter/gather functionality
    "LowPassEUR363Nov" -> new Target("EUR.nov2010", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, 99.0, lowPass, !exome, 363)
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
        if (!noIndels) add(new indelCall(target), new indelFilter(target), new indelEvaluation(target))
        add(new snpCall(target))
        add(new VQSR(target, !goldStandard))
        add(new applyVQSR(target, !goldStandard))
        add(new snpEvaluation(target))
      }
      if ( runGoldStandard ) {
        add(new VQSR(target, goldStandard))
        add(new applyVQSR(target, goldStandard))
      }
    }
  }


  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 4;
    phone_home = if ( LOCAL_ET ) GATKRunReport.PhoneHomeOption.STANDARD else GATKRunReport.PhoneHomeOption.AWS_S3
  }

  def bai(bam: File) = new File(bam + ".bai")

  // 1.) Unified Genotyper Base
  class GenotyperBase (t: Target) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 140
    this.nt = 2
    this.dcov = if ( t.isLowpass ) { 50 } else { 250 }
    this.stand_call_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.stand_emit_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.input_file :+= t.bamList
    this.D = new File(t.dbsnpFile)
  }

  // 1a.) Call SNPs with UG
  class snpCall (t: Target) extends GenotyperBase(t) {
    if (minimumBaseQuality >= 0)
      this.min_base_quality_score = minimumBaseQuality
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = t.rawVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (noBAQ ||  t.isExome) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY}
    this.analysisName = t.name + "_UGs"
    this.jobName =  queueLogDir + t.name + ".snpcall"
  }

  // 1b.) Call Indels with UG
  class indelCall (t: Target) extends GenotyperBase(t) {
    this.memoryLimit = 6
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
    this.V = t.rawIndelVCF
    this.out = t.filteredIndelVCF
    this.filterName ++= List("IndelQD", "IndelReadPosRankSum", "IndelFS")
    this.filterExpression ++= List("\"QD < 2.0\"", "\"ReadPosRankSum < -20.0\"", "\"FS > 200.0\"")
    if (t.nSamples >= 10) {
        this.filterName ++= List("IndelInbreedingCoeff")
        this.filterExpression ++= List("\"InbreedingCoeff < -0.8\"")
    }
    this.analysisName = t.name + "_VF"
    this.jobName =  queueLogDir + t.name + ".indelfilter"
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(t: Target, goldStandard: Boolean) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.nt = 2
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.input :+= ( if ( goldStandard ) { t.goldStandard_VCF } else { t.rawVCF } )
    this.resource :+= new TaggedFile( t.hapmapFile, "training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( omni_b37, "training=true,truth=true,prior=12.0" )
    this.resource :+= new TaggedFile( training_1000G, "training=true,prior=10.0" )
    this.resource :+= new TaggedFile( t.dbsnpFile, "known=true,prior=2.0" )
    this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS")
    if(t.nSamples >= 10) { // InbreedingCoeff is a population-wide statistic that requires at least 10 samples to calculate
        this.use_annotation ++= List("InbreedingCoeff")
    }
    if(!t.isExome) {
        this.use_annotation ++= List("DP")
    } else { // exome specific parameters
        this.resource :+= new TaggedFile( badSites_1000G, "bad=true,prior=2.0" )
        this.mG = 6
        if(t.nSamples <= 3) { // very few exome samples means very few variants
            this.mG = 4
            this.percentBad = 0.04
        }
    }
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile }
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.rscript_file = t.vqsrRscript
    this.analysisName = t.name + "_VQSR"
    this.jobName = queueLogDir + t.name + ".VQSR"
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  class applyVQSR (t: Target, goldStandard: Boolean) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 6
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.input :+= ( if ( goldStandard ) { t.goldStandard_VCF } else { t.rawVCF } )
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile}
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.ts_filter_level = t.trancheTarget
    this.out = t.recalibratedVCF
    this.analysisName = t.name + "_AVQSR"
    this.jobName = queueLogDir + t.name + ".applyVQSR"
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class EvalBase(t: Target) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.reference_sequence = t.reference
    this.comp :+= new TaggedFile(t.hapmapFile, "hapmap" )
    this.intervalsString ++= List(t.intervals)
    this.D = new File(t.dbsnpFile)
    this.sample = samples
  }

  // 5a.) SNP Evaluation (OPTIONAL) based on the cut vcf
  class snpEvaluation(t: Target) extends EvalBase(t) {
    if (t.reference == b37 || t.reference == hg19) this.comp :+= new TaggedFile( omni_b37, "omni" )
    this.eval :+= t.recalibratedVCF
    this.out =  t.evalFile
    this.analysisName = t.name + "_VEs"
    this.jobName = queueLogDir + t.name + ".snp.eval"
  }

  // 5b.) Indel Evaluation (OPTIONAL)
  class indelEvaluation(t: Target) extends EvalBase(t) {
    this.eval :+= t.filteredIndelVCF
    this.evalModule :+= "IndelStatistics"
    this.out =  t.evalIndelFile
    this.analysisName = t.name + "_VEi"
    this.jobName = queueLogDir + queueLogDir + t.name + ".indel.eval"
  }
}
