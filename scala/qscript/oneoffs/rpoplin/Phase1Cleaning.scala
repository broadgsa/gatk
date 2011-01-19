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

class Phase1Cleaning extends QScript {
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
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")
  private val dindelPilotCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
  private val dindelAFRCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/AFR.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelASNCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/ASN.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelEURCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/EUR.dindel_august_release.20110110.sites.vcf.gz"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","JPT","LWK","MXL","PUR","TSI","YRI")
  //val populations = List("ZZZ") // small set used for debugging
  
  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(2)
    this.jobTempDir = qscript.tmpDir
  }

  def script = {
    callThisChunk() // using scatter/gather capabilities of Queue so no need to for loop over 1Mb chunks of the chromosome
  }


  def callThisChunk() = {

    val interval = "%d".format(qscript.chr)
    for( population <- qscript.populations ) {
      val baseTmpName: String = qscript.outputTmpDir + "/" + population + ".phase1.chr" + qscript.chr.toString + "."
      val bamList: File = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.bam.lists/%s.chr%d.list".format(population, qscript.chr))
      val targetIntervals: File = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams/intervals/%s.chr%d.intervals".format(population, qscript.chr))

      // 1.) Create cleaning targets
      var target = new RealignerTargetCreator with CommandLineGATKArgs
      target.memoryLimit = Some(4)
      target.input_file :+= bamList
      target.intervalsString :+= interval
      target.out = targetIntervals
      target.mismatchFraction = Some(0.0)
      target.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
      target.rodBind :+= RodBind("indels1", "VCF", qscript.dindelPilotCalls)
      target.rodBind :+= RodBind("indels2", "VCF", qscript.dindelAFRCalls)
      target.rodBind :+= RodBind("indels3", "VCF", qscript.dindelEURCalls)
      target.rodBind :+= RodBind("indels4", "VCF", qscript.dindelASNCalls)
      target.jobName = baseName + population + ".target"

      // 2.) Clean without SW
      var clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + "cleaned.bam")
      clean.memoryLimit = Some(4)
      clean.input_file :+= bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = targetIntervals
      clean.out = cleanedBam
      clean.doNotUseSW = true
      clean.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE)
      clean.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
      clean.rodBind :+= RodBind("indels1", "VCF", qscript.dindelPilotCalls)
      clean.rodBind :+= RodBind("indels2", "VCF", qscript.dindelAFRCalls)
      clean.rodBind :+= RodBind("indels3", "VCF", qscript.dindelEURCalls)
      clean.rodBind :+= RodBind("indels4", "VCF", qscript.dindelASNCalls)
      clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true
      clean.jobName = baseName + population + ".clean"

      add(target, clean)
    }

  }
}