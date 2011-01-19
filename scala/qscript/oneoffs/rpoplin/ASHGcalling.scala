import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class ASHGcalling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="the chromosome to process", shortName="chr", required=true)
  var chr: Int = _

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/calls/"

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = "ALL.august"

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/humgen/gsa-hpprojects/august_cleaned_bams"

  private val tmpDir: File = new File("/broad/shptmp/rpoplin/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")
  private val targetIntervals: File = new File("/humgen/1kg/processing/allPopulations_wholeGenome_august_release/knownIndels.intervals")
  private val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
  private val dindelMask: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/pilot1.dindel.mask.bed"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("YRI","LWK","ASW","PUR","CEU","TSI","GBR","FIN","MXL","CHB","CHS","JPT")

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(2)
    this.DBSNP = qscript.dbSNP
    this.jobTempDir = qscript.tmpDir
  }

  class SamtoolsBaqFunction extends CommandLineFunction {
    @Input(doc="The input BAM file") var in_bam: File = _
    @Output(doc="The output BAM file") var out_bam: File = _
    def commandLine = "/humgen/gsa-scr1/rpoplin/samtools/samtools calmd -br %s %s > %s".format(in_bam.getAbsolutePath, qscript.reference, out_bam.getAbsolutePath)
  }

  class DeleteMeFunction extends CommandLineFunction {
    @Input(doc="The file to be deleted") var me: File = _
    @Input(doc="The file which must exist before we are allowed to delete") var trigger: File = _
    def commandLine = "rm -f %s".format(me.getAbsolutePath)
  }

  class DeleteMeAllFunction extends CommandLineFunction {
    @Input(doc="The file to be deleted") var me: File = _
    @Input(doc="The file which must exist before we are allowed to delete") var trigger: File = _
    def commandLine = "rm -f %s*".format(me.getAbsolutePath)
  }

  def script = {
    val basesPerJob: Int = 3000000
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

    /* CombineVariants parses the 800+ genotypes per record and is way too slow. Combine the vcf files together using grep, cat, and sortByRef.pl outside of Queue
    combineVariants = new CombineVariants with CommandLineGATKArgs
    combineVariants.rodBind = vcfChunks
    combineVariants.out = new TaggedFile(qscript.baseName + ".chr" + qscript.chr.toString + ".filtered.vcf", "vcf")
    combineVariants.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    combineVariants.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNSORTED)
    combineVariants.setKey = "null"
    add(combineVariants)
    */
  }

  def callThisChunk(interval: String, jobNumber: Int) = {

    val baseName: String = qscript.outputDir + "/chr" + qscript.chr.toString + "/" + qscript.baseName + ".chr" + qscript.chr.toString + "." + jobNumber.toString +"."
    var call = new UnifiedGenotyperV2 with CommandLineGATKArgs
    val rawCalls = new File(baseName + "raw.vcf")

    for( population <- qscript.populations ) {
      val baseTmpName: String = qscript.outputTmpDir + "/chr" + qscript.chr.toString + "/" + population + ".august.chr" + qscript.chr.toString + "." + jobNumber.toString +"."
      val bamList: File = new File("/humgen/1kg/processing/allPopulations_wholeGenome_august_release/bamLists/%s.chr%d.bam.list".format(population, qscript.chr))

      // 1.) Clean at known indels
      var clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + "cleaned.bam")
      clean.memoryLimit = Some(4)
      clean.input_file :+= bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = qscript.targetIntervals
      clean.out = cleanedBam
      clean.rodBind :+= RodBind("indels", "VCF", qscript.dindelCalls)
      clean.knownsOnly = true
      clean.LOD = Some(1.0)
      clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true
      clean.compress = Some(2)
      clean.jobName = baseName + population + ".clean"
      //clean.stripBam = true
      //clean.fileSystemUsage = "indium"

      // 2.) Apply BAQ calculation
      var baq = new SamtoolsBaqFunction
      val baqedBam = new File(baseTmpName + "cleaned.baq.bam")
      baq.memoryLimit = Some(4)
      baq.in_bam = cleanedBam
      baq.out_bam = baqedBam
      baq.jobName = baseName + population + ".baq"
      //baq.fileSystemUsage = "iodine"

      // 3a.) Delete cleaned bam
      var deleteClean = new DeleteMeFunction
      deleteClean.me = cleanedBam
      deleteClean.trigger = baqedBam
      deleteClean.jobName = baseName + population + ".deleteClean"
      //deleteClean.fileSystemUsage = "iodine"

      // 3b.) Index BAQ'ed bam
      var index = new SamtoolsIndexFunction
      index.bamFile = baqedBam
      index.jobName = baseName + population + ".index"
      //index.fileSystemUsage = "iodine"

      // 5a.) Delete BAQ'ed bam and index
      //var deleteBaq = new DeleteMeAllFunction
      //deleteBaq.me = baqedBam
      //deleteBaq.trigger = rawCalls
      //deleteBaq.jobName = baseName + population + ".deleteBaq"
      //deleteBaq.fileSystemUsage = "iodine"

      call.input_file :+= baqedBam

      //add(clean, baq, deleteClean, index, deleteBaq)
      add(clean, baq, deleteClean, index)
    }

    // 4.) Call with UGv2
    call.memoryLimit = Some(4)
    call.intervalsString :+= interval
    call.out = rawCalls
    call.dcov = Some(50)
    call.standard_min_confidence_threshold_for_calling = Some(50)
    call.standard_min_confidence_threshold_for_emitting = Some(30)
    call.min_mapping_quality_score = Some(20)
    call.min_base_quality_score = Some(20)
    call.pnrm = Some(org.broadinstitute.sting.playground.gatk.walkers.genotyper.AlleleFrequencyCalculationModel.Model.GRID_SEARCH)
    call.jobName = baseName + "call"
    //call.fileSystemUsage = "iodine"

    // 5b.) Filter near indels and HARD_TO_VALIDATE
    var filter = new VariantFiltration with CommandLineGATKArgs
    val filteredCalls = new File(baseName + "filtered.vcf")
    filter.memoryLimit = Some(1)
    filter.out = filteredCalls
    filter.intervalsString :+= interval
    filter.variantVCF = rawCalls
    filter.rodBind :+= RodBind("mask", "Bed", qscript.dindelMask)
    filter.maskName = "InDel"
    filter.filterName ++= List("HARD_TO_VALIDATE")
    filter.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    filter.jobName = baseName + "filter"
    //filter.fileSystemUsage = "indium"

    // 6.) Delete raw calls and index
    var deleteRawCalls = new DeleteMeAllFunction
    deleteRawCalls.me = rawCalls
    deleteRawCalls.trigger = filteredCalls
    deleteRawCalls.jobName = baseName + "deleteRawCalls"
    //deleteRawCalls.fileSystemUsage = "indium"

    add(call, filter, deleteRawCalls)
  }
}