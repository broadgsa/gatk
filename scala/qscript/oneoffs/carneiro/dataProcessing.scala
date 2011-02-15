
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.QScript


class dataProcessing extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="path to MarkDuplicates jar", shortName="dedup", required=true)
  var dedupJar: File = _

  @Input(doc="input BAM file", shortName="input", required=true)
  var inputBam: File = _

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = ""

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: File = _

  @Input(doc="the -L interval string to be used by GATK", shortName="L", required=false)
  var intervalString: String = ""

  @Input(doc="provide a .intervals file with the list of target intervals", shortName="intervals", required=false)
  var intervals: File = _


  def script = {

    val baseName: String        = swapExt(qscript.inputBam, ".bam", "").toString()
    def cleanedBam: String      = baseName + ".cleaned.bam"
    def dedupedBam: String      = baseName + ".cleaned.dedup.bam"
    def metricsFile: String     = swapExt(qscript.inputBam, "bam", "metrics").toString()
    def targetIntervals: String = baseName + ".indel.intervals"

    val reference: File       = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
    val dbSNP: File           = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")
    val dindelPilotCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
    val dindelAFRCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110126_dindel_august/AFR.dindel_august_release_merged_pilot1.20110126.sites.vcf.gz"
    val dindelASNCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110126_dindel_august/ASN.dindel_august_release_merged_pilot1.20110126.sites.vcf.gz"
    val dindelEURCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110126_dindel_august/EUR.dindel_august_release_merged_pilot1.20110126.sites.vcf.gz"

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.jarFile = qscript.gatkJar
      this.reference_sequence = reference
      this.memoryLimit = Some(2)
      this.jobTempDir = qscript.outputTmpDir
    }


    val target = new RealignerTargetCreator with CommandLineGATKArgs
    target.memoryLimit = Some(4)
    target.input_file :+= qscript.inputBam
    target.out = new File(targetIntervals)
    target.mismatchFraction = Some(0.0)
    target.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    target.rodBind :+= RodBind("indels1", "VCF", dindelPilotCalls)
    target.rodBind :+= RodBind("indels2", "VCF", dindelAFRCalls)
    target.rodBind :+= RodBind("indels3", "VCF", dindelEURCalls)
    target.rodBind :+= RodBind("indels4", "VCF", dindelASNCalls)
    target.jobName = baseName + ".target"
    if (!qscript.intervalString.isEmpty()) target.intervalsString ++= List(qscript.intervalString)
    if (qscript.intervals != Nil) target.intervals ++= List(qscript.intervals)

    // 2.) Clean without SW

    val clean = new IndelRealigner with CommandLineGATKArgs
    clean.memoryLimit = Some(4)
    clean.input_file :+= qscript.inputBam
    clean.targetIntervals = new File(targetIntervals)
    clean.out = new File(cleanedBam)
    clean.doNotUseSW = true
    clean.baq = Some(org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE)
    clean.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    clean.rodBind :+= RodBind("indels1", "VCF", dindelPilotCalls)
    clean.rodBind :+= RodBind("indels2", "VCF", dindelAFRCalls)
    clean.rodBind :+= RodBind("indels3", "VCF", dindelEURCalls)
    clean.rodBind :+= RodBind("indels4", "VCF", dindelASNCalls)
    clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true
    clean.jobName = baseName + ".clean"
    if (!qscript.intervalString.isEmpty()) clean.intervalsString ++= List(qscript.intervalString)
    if (qscript.intervals != Nil) clean.intervals ++= List(qscript.intervals)

    // 3.) Mark Duplicates
    val dedup = new PicardBamJarFunction{
      @Input(doc="cleaned bam") var clean: File = new File(cleanedBam)
      @Output(doc="deduped bam") var dedup: File = new File(dedupedBam)
      override def inputBams = List(clean)
      override def outputBam = dedup
      override def commandLine = super.commandLine + " M=" + metricsFile
      sortOrder = null
    }
    dedup.memoryLimit = Some(8)
    dedup.jarFile = qscript.dedupJar
    dedup.jobName = baseName + ".dedup"

  //  val cov = new CountCovariates with CommandLineGATKArgs


    add(target, clean, dedup)
  }
}