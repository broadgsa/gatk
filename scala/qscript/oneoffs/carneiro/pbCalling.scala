import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class pbCalling extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=true)
  var outputDir: String = "./"


  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil


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
          val isLowpass: Boolean,
          val isCCS: Boolean) {
    val name = qscript.outputDir + baseName
    val clusterFile = new File(name + ".clusters")
    val rawVCF = new File(name + ".raw.vcf")
    val filteredVCF = new File(name + ".filtered.vcf")
    val titvRecalibratedVCF = new File(name + ".titv.recalibrated.vcf")
    val titvTranchesFile = new File(name + ".titv.tranches")
    val recalibratedVCF = new File(name + ".ts.recalibrated.vcf")
    val tranchesFile = new File(name + ".ts.tranches")
    val cutVCF = new File(name + ".cut.vcf")
    val evalFile = new File(name + ".eval")
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
  val dbSNP_b37_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.rod"              // Special case for NA12878 collections that can't use 132 because they are part of it.
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf"
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
  val ccs: Boolean = true

  val targetDataSets: Map[String, Target] = Map(
    "HiSeq" -> new Target("NA12878.HiSeq", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals", 2.07, !lowPass, !ccs),
    "FIN" ->  new Target("FIN", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass, !ccs),
    "WEx" ->  new Target("NA12878.WEx", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, !lowPass, !ccs),
    "TGPWExGdA" -> new Target("1000G.WEx.GdA", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"),        // BUGBUG: reduce from 60 to 20 people
              new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass, !ccs),
    "LowPassN60" -> new Target("lowpass.N60", b36, dbSNP_b36, hapmap_b36, indelMask_b36,
              new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
              new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, lowPass, !ccs),          // chunked interval list to use with Queue's scatter/gather functionality
    "LowPassAugust" -> new Target("ALL.august.v4", b37, dbSNP_b37, hapmap_b37, indelMask_b37,                               // BUGBUG: kill this, it is too large
              new File("/humgen/1kg/processing/allPopulations_chr20_august_release.cleaned.merged.bams/ALL.cleaned.merged.list"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass, !ccs),
    "LowPassEUR363Nov" -> new Target("EUR.nov2010", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass, !ccs),
    "WExTrio" -> new Target("NA12878Trio.WEx", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.bwa.cleaned.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/trio/snps/NA12878Trio.WEx.filtered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass, !ccs),
    "pacbio" -> new Target("pacbio", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/data/pacbio.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/analisys/snps/pacbio.filtered.vcf"),
				      "/humgen/gsa-scr1/carneiro/prj/pacbio/data/pacbio.hg19.intervals", 1.8, !lowPass, !ccs),
    "pb200" -> new Target("pb200", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/data/pb200.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/analisys/snps/pb200.filtered.vcf"),
				      "/humgen/gsa-scr1/carneiro/prj/pacbio/data/pb200.hg19.intervals", 1.8, !lowPass, !ccs),
    "pb2k" -> new Target("pb2k", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/data/pb2k.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/analisys/snps/pb2k.filtered.vcf"),
				      "/humgen/gsa-scr1/carneiro/prj/pacbio/data/pb2k.hg19.intervals", 1.8, !lowPass, !ccs),
    "cc200" -> new Target("cc200", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/data/cc200.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/analisys/snps/cc200.filtered.vcf"),
				      "/humgen/gsa-scr1/carneiro/prj/pacbio/data/cc200.hg19.intervals", 1.8, !lowPass, ccs),
    "cc2k" -> new Target("cc2k", b37, dbSNP_b37_129, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/data/cc2k.recal.bam"),
              new File("/humgen/gsa-scr1/carneiro/prj/pacbio/analisys/snps/cc2k.filtered.vcf"),
				      "/humgen/gsa-scr1/carneiro/prj/pacbio/data/cc2k.hg19.intervals", 1.8, !lowPass, ccs)
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
      add(new UnifiedGenotyper(target))
      add(new VariantFiltration(target))
      add(new VQSR(target, !goldStandard))
      add(new applyVQSR(target, !goldStandard))
      add(new VariantCut(target))
      add(new VariantEvaluation(target))
    }
  }

  def bai(bam: File) = new File(bam + ".bai")

  val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

  // 1.) Call SNPs with UG
  class UnifiedGenotyper(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper {
    this.jarFile = gatkJarFile
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
    this.dcov = if ( t.isLowpass ) { 50 } else { 250 }
    this.stand_call_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.stand_emit_conf = if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.input_file :+= t.bamList
    this.out = t.rawVCF
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.analysisName = t.name + "_UG"
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
    // Ridiculous workaround to get pacbio data to run.. never commit this!
    this.deletions = 0.5
    this.mbq = 10
  }

  // 2.) Filter SNPs
  class VariantFiltration(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantFiltration {
    this.jarFile = gatkJarFile
    this.reference_sequence = t.reference
    this.intervalsString ++= List(t.intervals)
    this.scatterCount = 10
    this.variantVCF = t.rawVCF
    this.out = t.filteredVCF
    this.filterName ++= List("HARD_TO_VALIDATE")
    this.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    this.analysisName = t.name + "_VF"
  }

  class VQSR(t: Target, goldStandard: Boolean) extends VariantRecalibrator {
    this.memoryLimit = 6
    this.intervalsString ++= List(t.intervals)
    this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
    this.rodBind :+= RodBind("hapmap", "VCF", t.hapmapFile)
    if( t.hapmapFile.contains("b37") )
      this.rodBind :+= RodBind("1kg", "VCF", omni_b37)
    else if( t.hapmapFile.contains("b36") )
      this.rodBind :+= RodBind("1kg", "VCF", omni_b36)
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
    this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile }
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.allPoly = true
    this.tranche ++= List("0.1", "0.5", "0.7", "1.0", "1.1", "1.2", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "2.1", "2.2", "2.5","3.0", "5.0", "10.0")
  }

  class applyVQSR (t: Target, goldStandard: Boolean) extends ApplyRecalibration {
    this.memoryLimit = 4
    this.intervalsString ++= List(t.intervals)
    this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
    this.tranches_file = if ( goldStandard ) { t.goldStandardTranchesFile } else { t.tranchesFile}
    this.recal_file = if ( goldStandard ) { t.goldStandardRecalFile } else { t.recalFile }
    this.fdr_filter_level = 2.0
    this.out = t.recalibratedVCF
  }

  // 5.) Variant Cut filter out the variants marked by recalibration to the 99% tranche
  class VariantCut(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.ApplyVariantCuts {
      this.jarFile = gatkJarFile
      this.reference_sequence = t.reference
      this.rodBind :+= RodBind("input", "VCF",  t.recalibratedVCF )
      this.analysisName = t.name + "_VC"
      this.intervalsString ++= List(t.intervals)
      this.out = t.cutVCF
      this.tranchesFile = t.tranchesFile
      this.fdr_filter_level = 1.0
      if (t.dbsnpFile.endsWith(".rod"))
        this.DBSNP = new File(t.dbsnpFile)
	    else if (t.dbsnpFile.endsWith(".vcf"))
        this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }

  // 6.) Variant Evaluation  based on the sensitivity recalibrated vcf
  class VariantEvaluation(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantEval {
      this.jarFile = gatkJarFile
      val name: String = t.name
      this.reference_sequence = t.reference
      this.rodBind :+= RodBind("comp", "VCF", t.hapmapFile)
      this.rodBind :+= RodBind("eval", "VCF", t.cutVCF)
      this.analysisName = name + "_VE"
      this.intervalsString ++= List(t.intervals)
      this.EV ++= List("GenotypeConcordance")
      this.out =  t.evalFile
    // Ridiculous workaround to get pacbio data to run.. never commit this!
      this.sample ++= List("NA12878")
  }
}
