import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class Phase1Calling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="the chromosome to process", shortName="chr", required=false)
  var chr: Int = 20

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = "/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams"

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams"

  private val tmpDir: File = new File("/broad/shptmp/rpoplin/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")
  private val dindelPilotCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
  private val dindelAFRCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/AFR.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelASNCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/ASN.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelEURCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/EUR.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelMask: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/pilot1.dindel.mask.bed"
  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val g1k = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf"
  val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","JPT","LWK","MXL","PUR","TSI","YRI")    
  //val populations = List("JPT","ASN","AMR")
  //val populations = List("EUR","AMR","ASN","AFR")
  //val populations = List("FIN", "LWK")
  private val intervals: String = "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals"
  //val populations = List("ZZZ") // small set used for debugging

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(3)
    this.jobTempDir = qscript.tmpDir
    this.DBSNP = qscript.dbSNP
  }

  def script = {
    callThisChunk() // using scatter/gather capabilities of Queue so no need to for loop over 1Mb chunks of the chromosome
  }

  def callThisChunk() = {

    val interval = "%d".format(qscript.chr)
    for( population <- qscript.populations ) {
      val baseName: String = qscript.outputDir + "/" + population + ".phase1.chr" + qscript.chr.toString
      var bamList: File = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams/%s.phase1.chr%d.cleaned.bam".format(population, qscript.chr))
      if( population == "ASN" || population == "EUR" || population == "AFR" || population == "AMR" ) {
        bamList = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams/%s.chr%d.cleaned.list".format(population, qscript.chr))
      }

      val rawCalls = new File(baseName + ".raw.vcf")
      val filteredCalls = new File(baseName + ".filtered.vcf")
      val clusterFile = new File(baseName + ".omni.clusters")
      val recalibratedCalls = new File(baseName + ".recal.vcf")
      val tranchesFile = new File(baseName + ".ts.omni.tranches")

      var call = new UnifiedGenotyper with CommandLineGATKArgs
      call.intervalsString ++= List(qscript.intervals)
      call.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
      call.dcov = Some( 50 )
      call.stand_call_conf = Some( 4.0 )
      call.stand_emit_conf = Some( 4.0 )
      call.input_file :+= bamList
      call.out = rawCalls
      call.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY)
      call.analysisName = baseName + "_UG"

      var filter = new VariantFiltration with CommandLineGATKArgs
      filter.intervalsString ++= List(qscript.intervals)
      filter.scatterCount = 10
      filter.variantVCF = rawCalls
      filter.out = filteredCalls
      filter.filterName ++= List("HARD_TO_VALIDATE")
      filter.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
      filter.analysisName = baseName + "_VF"
      //filter.rodBind :+= RodBind("mask", "Bed", qscript.dindelMask)
      //filter.maskName = "InDel"

      var gvc = new GenerateVariantClusters with CommandLineGATKArgs
      gvc.rodBind :+= RodBind("hapmap", "VCF", qscript.hapmap)
      gvc.rodBind :+= RodBind("1kg", "VCF", qscript.omni)
      gvc.rodBind :+= RodBind("input", "VCF", filteredCalls )
      gvc.clusterFile = clusterFile
      gvc.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      gvc.analysisName = baseName + "_GVC"
      gvc.intervalsString ++= List(qscript.intervals)
      gvc.qual = Some(100) // clustering parameters to be updated soon pending new experimentation results
      gvc.std = Some(4.5)
      gvc.mG = Some(6)

      var vr = new VariantRecalibrator with CommandLineGATKArgs
      vr.rodBind :+= RodBind("1kg", "VCF", qscript.omni)
      vr.rodBind :+= RodBind("hapmap", "VCF", qscript.hapmap)
      vr.rodBind :+= RodBind("truthOmni", "VCF", qscript.omni)
      vr.rodBind :+= RodBind("truthHapMap", "VCF", qscript.hapmap)
      vr.rodBind :+= RodBind("input", "VCF", filteredCalls )
      vr.clusterFile = clusterFile
      vr.analysisName = baseName + "_VR"
      vr.intervalsString ++= List(qscript.intervals)
      vr.ignoreFilter ++= List("HARD_TO_VALIDATE")
      vr.target_titv = Some(2.3)
      vr.sm = Some(org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY)
      vr.tranche ++= List("0.1", "1.0", "2.0", "3.0", "5.0", "10.0", "100.0")
      vr.out = recalibratedCalls
      vr.priorDBSNP = Some(10.0)
      vr.priorHapMap = Some(12.0)
      vr.prior1KG = Some(12.0)
      vr.tranchesFile = tranchesFile      

      add(call, filter, gvc, vr)
    }

  }
}