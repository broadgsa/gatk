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

class Phase1WholeGenome extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="the chromosome to process", shortName="chr", required=true)
  var chr: Int = _

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=true)
  var outputTmpDir: String = _

  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")
  private val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS","JPT","LWK","MXL","PUR","TSI","YRI") //,"PPN")

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(3)
    this.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )
  }

  class AnalysisPanel(val baseName: String, val pops: List[String], val jobNumber: Int) {
    val rawVCFsnps = new File(qscript.outputDir + "/calls/chr" + qscript.chr.toString + "/" + baseName + "/" + baseName + ".phase1.chr" + qscript.chr.toString + "." + jobNumber + ".raw.snps.vcf")
    val rawVCFindels = new File(qscript.outputDir + "/calls/chr" + qscript.chr.toString + "/" + baseName + "/" + baseName + ".phase1.chr" + qscript.chr.toString + "." + jobNumber + ".raw.indels.vcf")

    val callSnps = new UnifiedGenotyper with CommandLineGATKArgs
    callSnps.out = rawVCFsnps
    callSnps.dcov = Some( 50 )
    callSnps.stand_call_conf = Some( 4.0 )
    callSnps.stand_emit_conf = Some( 4.0 )
    callSnps.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY)
    callSnps.jobName = qscript.outputTmpDir + "/calls/chr" + qscript.chr.toString + "/" +baseName + ".phase1.chr" + qscript.chr.toString + "." + jobNumber + ".raw.snps"
    callSnps.exactCalculation = Some(org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel.ExactCalculation.LINEAR_EXPERIMENTAL)
    
    val callIndels = new UnifiedGenotyper with CommandLineGATKArgs
    callIndels.out = rawVCFindels
    callIndels.dcov = Some( 50 )
    callIndels.stand_call_conf = Some( 10.0 )
    callIndels.stand_emit_conf = Some( 10.0 )
    callIndels.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF)
    callIndels.glm = Some(org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.DINDEL)
    callIndels.minIndelCnt = Some(5)
    callIndels.read_filter :+= "Platform454"
    callIndels.jobName = qscript.outputTmpDir + "/calls/chr" + qscript.chr.toString + "/" +baseName + ".phase1.chr" + qscript.chr.toString + "." + jobNumber + ".raw.indels"
    callIndels.exactCalculation = Some(org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel.ExactCalculation.LINEAR_EXPERIMENTAL)
    callIndels.abort_at_too_much_coverage = Some(4500)
  }

  def script = {
    val basesPerJob: Int = 700000
    val lastBase: Int = qscript.chromosomeLength(qscript.chr - 1)
    var start: Int = 1
    var stop: Int = start - 1 + basesPerJob
    if( stop > lastBase ) { stop = lastBase }
    var jobNumber: Int = 1
    while( jobNumber < (lastBase.toFloat / basesPerJob.toFloat) + 1.0) {
      callThisChunk("%d:%d-%d".format(qscript.chr, start, stop), jobNumber)
      start += basesPerJob
      stop += basesPerJob
      if( stop > lastBase ) { stop = lastBase }
      jobNumber += 1
    }
  }


  def callThisChunk(interval: String, jobNumber: Int) = {

    val AFR = new AnalysisPanel("AFR", List("LWK","YRI","ASW","CLM","PUR"), jobNumber)
    val AMR = new AnalysisPanel("AMR", List("MXL","CLM","PUR","ASW"), jobNumber)
    val EUR = new AnalysisPanel("EUR", List("CEU","FIN","GBR","TSI","IBS","MXL","CLM","PUR","ASW"), jobNumber)
    val ASN = new AnalysisPanel("ASN", List("CHB","CHS","JPT","MXL","CLM","PUR"), jobNumber)
    //val PAA = new AnalysisPanel("PAA", List("PPN","YRI","CHB"), jobNumber)
    val analysisPanels = List(AFR, ASN, AMR, EUR)

    for( population <- qscript.populations ) {
      val baseTmpName: String = qscript.outputTmpDir + "/calls/chr" + qscript.chr.toString + "/" + population + ".phase1.chr" + qscript.chr.toString + "." + jobNumber.toString + "."
      val bamList: File = new File("/humgen/1kg/processing/allPopulations_wholeGenome_phase1_release/bam_lists/%s.list".format(population))
      val targetIntervals: File = new File(baseTmpName + "target.intervals")

      // 1.) Create cleaning targets
      val target = new RealignerTargetCreator with CommandLineGATKArgs
      target.memoryLimit = Some(3)
      target.input_file :+= bamList
      target.intervalsString :+= interval
      target.out = targetIntervals
      target.mismatchFraction = Some(0.0)
      target.maxIntervalSize = Some(700)
      target.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      target.jobName = baseTmpName + "target"
      target.isIntermediate = true

      // 2.) Clean without SW
      val clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + "cleaned.bam")
      clean.memoryLimit = Some(4)
      clean.input_file :+= bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = targetIntervals
      clean.out = cleanedBam
      clean.doNotUseSW = true
      clean.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE)
      clean.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      clean.jobName = baseTmpName + "clean"
      clean.isIntermediate = true
      clean.compress = Some(0)
      clean.index_output_bam_on_the_fly = Some(true)
      clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true

      add(target, clean)

      for( a <- analysisPanels ) {
        for( p <- a.pops) {
          if( p == population ) {
            a.callSnps.input_file :+= cleanedBam
            a.callIndels.input_file :+= cleanedBam
          }
        }
      }
    }

    for( a <- analysisPanels ) {
      a.callSnps.intervalsString :+= interval
      a.callIndels.intervalsString :+= interval
      add(a.callSnps, a.callIndels)
    }

  }
}
