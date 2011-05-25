import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class Phase1ProjectConsensus extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Input(doc="the chromosome to process", shortName="onlyChr20", required=false)
  var onlyChr20: Boolean = false

  @Input(doc="the chromosome to process", shortName="indelsOnly", required=false)
  var indelsOnly: Boolean = false

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=true)
  var outputTmpDir: String = "/broad/shptmp/delangel"

  @Input(doc="Generate bam files", shortName="generateBAMs", required=false)
   var generateBAMs: Boolean = false

  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")
  private val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)
  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS","JPT","LWK","MXL","PUR","TSI","YRI")
  private val snpAlleles: String = "/humgen/1kg/processing/production_wgs_phase1/consensus/ALL.phase1.wgs.union.pass.sites.vcf"
  private val indelAlleles: String = "//humgen/1kg/processing/production_wgs_phase1/consensus/ALL.indels.combined.chr20.vcf"

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(3)
    this.jobQueue = "gsa"

  }

  class AnalysisPanel(val baseName: String, val pops: List[String], val jobNumber: Int, val chr: String) {
    val rawVCFsnps = new File(qscript.outputDir + "/calls/chr" + chr + "/" + baseName + "/" + baseName + ".phase1.chr" + chr + "." + jobNumber + ".raw.snps.vcf")
    val rawVCFindels = new File(qscript.outputDir + "/calls/chr" + chr + "/" + baseName + "/" + baseName + ".phase1.chr" + chr + "." + jobNumber + ".raw.indels.vcf")

    val callSnps = new UnifiedGenotyper with CommandLineGATKArgs
    callSnps.out = rawVCFsnps
    callSnps.dcov = 50
    callSnps.stand_call_conf = 4.0
    callSnps.stand_emit_conf = 4.0
    callSnps.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
    callSnps.jobName = qscript.outputTmpDir + "/calls/chr" + chr + "/" +baseName + ".phase1.chr" + chr + "." + jobNumber + ".raw.snps"
    callSnps.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    callSnps.genotyping_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    //callSnps.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
    callSnps.rodBind :+= RodBind("alleles", "VCF", qscript.snpAlleles)
    callSnps.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )
    callSnps.sites_only = true

    val callIndels = new UnifiedGenotyper with CommandLineGATKArgs
    callIndels.out = rawVCFindels
    callIndels.dcov = 50
    callIndels.stand_call_conf = 4.0
    callIndels.stand_emit_conf = 4.0
    callIndels.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    callIndels.jobName = qscript.outputTmpDir + "/calls/chr" + chr + "/" +baseName + ".phase1.chr" + chr + "." + jobNumber + ".raw.indels"
    callIndels.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    callIndels.genotyping_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    callIndels.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
    callIndels.rodBind :+= RodBind("alleles", "VCF", qscript.indelAlleles)
    callIndels.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )
    callIndels.sites_only = true
    // indel specific stuff: add FisherStrand, add F-test statistic (excessive hets/homs)

    callIndels.A ++= List ("FisherStrand","GLstats")
  }

  class Chromosome(val inputChr: Int) {
    var chr: String = inputChr.toString
    if(inputChr == 23) { chr = "X" }

    val snpCombine = new CombineVariants with CommandLineGATKArgs
    val snpChrVCF = new File(qscript.outputDir + "/calls/" + "combined.phase1.chr" + chr + ".raw.snps.vcf")
    snpCombine.out = snpChrVCF
    snpCombine.intervalsString :+= chr
    val indelCombine = new CombineVariants with CommandLineGATKArgs
    val indelChrVCF = new File(qscript.outputDir + "/calls/" + "combined.phase1.chr" + chr + ".raw.indels.vcf")
    indelCombine.out = indelChrVCF
    indelCombine.intervalsString :+= chr
  }

  def script = {

    var chrList =  List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    if (qscript.onlyChr20) {
      chrList = List(20)
    }
    for(chr <- chrList) {
      val chrObject = new Chromosome(chr)
      val basesPerJob: Int = 1000000
      val lastBase: Int = qscript.chromosomeLength(chr - 1)
      var start: Int = 1
      var stop: Int = start - 1 + basesPerJob
      if( stop > lastBase ) { stop = lastBase }
      var jobNumber: Int = 1
      while( jobNumber < (lastBase.toFloat / basesPerJob.toFloat) + 1.0) {
        if( chr != 23 ) {
          callThisChunk("%d:%d-%d".format(chr, start, stop), jobNumber, chr, chrObject)
        } else {
          callThisChunk("X:%d-%d".format(start, stop), jobNumber, chr, chrObject)
        }
        start += basesPerJob
        stop += basesPerJob
        if( stop > lastBase ) { stop = lastBase }
        jobNumber += 1
      }
      add(chrObject.indelCombine)
      if (!qscript.indelsOnly)
        add(chrObject.snpCombine)
    }
  }


  def callThisChunk(interval: String, jobNumber: Int, inputChr: Int, chrObject: Chromosome) = {

    var chr: String = inputChr.toString
    if(inputChr == 23) { chr = "X" }

    val AFRadmix = new AnalysisPanel("AFR.admix", List("LWK","YRI","ASW","CLM","PUR"), jobNumber, chr)
    val AMRadmix = new AnalysisPanel("AMR.admix", List("MXL","CLM","PUR","ASW"), jobNumber, chr)
    val EURadmix = new AnalysisPanel("EUR.admix", List("CEU","FIN","GBR","TSI","IBS","MXL","CLM","PUR","ASW"), jobNumber, chr)
    val ASNadmix = new AnalysisPanel("ASN.admix", List("CHB","CHS","JPT","MXL","CLM","PUR"), jobNumber, chr)
    val AFR = new AnalysisPanel("AFR", List("LWK","YRI","ASW"), jobNumber, chr)
    val AMR = new AnalysisPanel("AMR", List("MXL","CLM","PUR"), jobNumber, chr)
    val EUR = new AnalysisPanel("EUR", List("CEU","FIN","GBR","TSI","IBS"), jobNumber, chr)
    val ASN = new AnalysisPanel("ASN", List("CHB","CHS","JPT"), jobNumber, chr)
    val ALL = new AnalysisPanel("ALL", List("LWK","YRI","ASW","MXL","CLM","PUR","CEU","FIN","GBR","TSI","IBS","CHB","CHS","JPT"), jobNumber, chr)

    val analysisPanels = List(AFR, ASN, AMR, EUR, AFRadmix, ASNadmix, AMRadmix, EURadmix) //ALL

    val snpCombine = new CombineVariants with CommandLineGATKArgs
    val indelCombine = new CombineVariants with CommandLineGATKArgs
    val combinedIndelChunk = new File(qscript.outputDir + "/calls/chr" + chr + "/" + "combined.phase1.chr" + chr + "." + jobNumber + ".raw.indels.vcf")
    val combinedSnpChunk = new File(qscript.outputDir + "/calls/chr" + chr + "/" + "combined.phase1.chr" + chr + "." + jobNumber + ".raw.snps.vcf")

    indelCombine.out = combinedIndelChunk
    indelCombine.jobName = qscript.outputTmpDir + "/calls/chr" + chr + "/" + "combined.phase1.chr" + chr + "." + jobNumber + ".raw.indels"
    indelCombine.intervalsString :+= interval
    indelCombine.mergeInfoWithMaxAC = true
    indelCombine.priority = "AFR.admix,AMR.admix,EUR.admix,ASN.admix,AFR,AMR,EUR,ASN" //ALL,

    snpCombine.out = combinedSnpChunk
    snpCombine.jobName = qscript.outputTmpDir + "/calls/chr" + chr + "/" + "combined.phase1.chr" + chr + "." + jobNumber + ".raw.snps"
    snpCombine.intervalsString :+= interval
    snpCombine.mergeInfoWithMaxAC = true
    snpCombine.priority = "AFR.admix,AMR.admix,EUR.admix,ASN.admix,AFR,AMR,EUR,ASN" //ALL,

    for( population <- qscript.populations ) {
      val baseTmpName: String = qscript.outputTmpDir + "/calls/chr" + chr + "/" + population + ".phase1.chr" + chr + "." + jobNumber.toString + "."
      val bamList: File = new File("/humgen/1kg/processing/production_wgs_phase1/bam_lists/%s.list".format(population))
      val targetIntervals: File = new File(baseTmpName + "target.intervals")

      // 1.) Create cleaning targets
      val target = new RealignerTargetCreator with CommandLineGATKArgs
      target.memoryLimit = 4
      target.input_file :+= bamList
      target.intervalsString :+= interval
      target.out = targetIntervals
      target.mismatchFraction = 0.0
      target.maxIntervalSize = 700
      target.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      target.jobName = baseTmpName + "target"
      //target.isIntermediate = true
      target.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )

      // 2.) Clean without SW
      val clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + "cleaned.bam")
      clean.memoryLimit = 6
      clean.input_file :+= bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = targetIntervals
      clean.out = cleanedBam
      clean.doNotUseSW = true
      clean.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
      clean.simplifyBAM = true
      clean.rodBind :+= RodBind("indels1", "VCF", qscript.dindelCalls)
      clean.jobName = baseTmpName + "clean"
      //clean.isIntermediate = true
      clean.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )

      if (qscript.generateBAMs) {
        add(target, clean)
      }

      for( a <- analysisPanels ) {
        for( p <- a.pops) {
          if( p == population ) {
            a.callIndels.input_file :+= cleanedBam
            a.callSnps.input_file :+= cleanedBam
          }
        }
      }
    }

    for( a <- analysisPanels ) {
      a.callSnps.intervalsString :+= interval
      a.callIndels.intervalsString :+= interval
      if(a.baseName == "ALL") {
        a.callIndels.memoryLimit = 4
        a.callSnps.memoryLimit = 4
      }
      if (!qscript.indelsOnly)
        add(a.callSnps)
      add(a.callIndels)

      snpCombine.rodBind :+= RodBind(a.baseName, "VCF", a.callSnps.out)
      indelCombine.rodBind :+= RodBind(a.baseName, "VCF", a.callIndels.out)
    }

    add(indelCombine)
    if (!qscript.indelsOnly)
      add(snpCombine)

    chrObject.snpCombine.rodBind :+= RodBind("ALL" + jobNumber.toString, "VCF", snpCombine.out)
    chrObject.indelCombine.rodBind :+= RodBind("ALL" + jobNumber.toString, "VCF", indelCombine.out)
  }
}
