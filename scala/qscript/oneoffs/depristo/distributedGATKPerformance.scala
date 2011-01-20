import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;

class DistributedGATKPerformance extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = ""

  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil

  @Argument(shortName="long", doc="runs long calculations", required=false)
  var long: Boolean = false


  //@Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

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
    def rawVCF(part: String) = new File(name + "." + part + ".raw.vcf")
    val filteredVCF = new File(name + ".filtered.vcf")
    val titvRecalibratedVCF = new File(name + ".titv.recalibrated.vcf")
    val tsRecalibratedVCF = new File(name + ".ts.recalibrated.vcf")
    val goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    val goldStandardClusterFile = new File(goldStandardName + ".clusters")
  }

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
              "chr1", 2.07, !lowPass),
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
        new File("/humgen/gsa-scr1/carneiro/prj/trio/NA12878Trio.WEx.hg19.bam"),
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

    val nWays = if (long) List(1, 2, 5, 10) else List(25, 50, 100)

    for (target <- targets) {
        for (nWaysParallel <- nWays) {
//        for (nWaysParallel <- List(2, 5)) {
            val aname = "distN" + nWaysParallel;
            val coordinationFile = new File(target.name + "." + aname + ".distributed.txt")
            for ( part <- 1 to nWaysParallel) {
                add(new UnifiedGenotyper(target, coordinationFile, part, aname + ".part" + part))
            }
        }
    }
  }

  // 1.) Call SNPs with UG
  class UnifiedGenotyper(t: Target, coordinationFile: File, i: Int, part: String) extends org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.processingTracker = coordinationFile
    this.intervalsString ++= List(t.intervals)
    this.dcov = Some( if ( t.isLowpass ) { 50 } else { 250 } )
    this.stand_call_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.stand_emit_conf = Some( if ( t.isLowpass ) { 4.0 } else { 30.0 } )
    this.input_file :+= t.bamList
    this.out = t.rawVCF(part)
    this.baq = Some( if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE})
    this.analysisName = t.name + "_UG." + part
    if ( i == 1 ) this.performanceLog = new File(coordinationFile.getAbsolutePath + "." + part + ".pf.log")
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }
}
