package oneoffs.depristo

//import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction


class CleaningTest extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=false)
  var gatkJar: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/GenomeAnalysisTK.jar")

  @Input(doc="the chromosome to process", shortName="chr", required=false)
  var chr: String = "20"

  @Input(doc="the chromosome to process", shortName="L", required=false)
  var range: String = _

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/isizeConstrainedRealigner/"

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/broad/shptmp/depristo/tmp"

  @Input(doc="path to Picard FixMateInformation.jar.  See http://picard.sourceforge.net/ .", required=false)
  var picardFixMatesJar: File = new java.io.File("/seq/software/picard/current/bin/FixMateInformation.jar")
  var picardValidateJar: File = new java.io.File("/seq/software/picard/current/bin/ValidateSamFile.jar")
  var picardSortSamJar: File = new java.io.File("/seq/software/picard/current/bin/SortSam.jar")

  private val tmpDir: File = new File("/broad/shptmp/depristo/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")
  private val dindelEURCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/EUR.dindel_august_release.20110110.sites.vcf.gz"
//  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)

//  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 4
    this.jobTempDir = qscript.tmpDir
  }

  def script = {
    val interval = qscript.chr
    val bamList: File = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/isizeConstrainedRealigner/CEU.chr%s.list".format(qscript.chr))
    //val bamList: File = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/isizeConstrainedRealigner/FIN.chr%s.3samples.list".format(qscript.chr))
    val targetIntervals: File = new File("%s/chr_%s.intervals".format(outputDir, qscript.chr))

    Console.println("interval " + interval)

    // 1.) Create cleaning targets
    var target = new RealignerTargetCreator with CommandLineGATKArgs
    target.input_file :+= bamList
    target.intervalsString :+= interval
    target.out = targetIntervals
    target.mismatchFraction = 0.0
    target.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
    target.rodBind :+= RodBind("indels3", "VCF", qscript.dindelEURCalls)
    //target.jobName = baseName + ".target"
    add(target)

    for ( cm <- List(true, false) ) {
      // 2.) Clean without SW
      var clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(outputDir + "cleaned.cm_%b.bam".format(cm))

      clean.input_file :+= bamList
      clean.intervalsString :+= interval + (if ( range != null ) ":" + range else "")
      clean.targetIntervals = targetIntervals
      clean.out = if ( cm ) cleanedBam else new File(cleanedBam + ".intermediate.bam")
      clean.doNotUseSW = true
      clean.constrainMovement = cm
      clean.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
      clean.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP)
      clean.rodBind :+= RodBind("indels3", "VCF", qscript.dindelEURCalls)
      //clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true
      //clean.jobName = baseName + cm + ".clean"

      Console.println("CLEAN")
      add(clean)

      if ( ! cm ) {
          // Explicitly run fix mates if the function won't be scattered.
          val fixMates = new PicardBamFunction {
            // Declare inputs/outputs for dependency tracking.
            @Input(doc="unfixed bam") var unfixed: File = _
            @Output(doc="fixed bam") var fixed: File = _
            def inputBams = List(unfixed)
            def outputBam = fixed
          }

          //fixMates.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
          fixMates.memoryLimit = 4
          fixMates.jarFile = qscript.picardFixMatesJar
          fixMates.unfixed = clean.out
          fixMates.fixed = cleanedBam
          //fixMates.analysisName = "FixMates"

          // Add the fix mates explicitly
          Console.println("fixMates")
          add(fixMates)
      }

      val validate = new JavaCommandLineFunction {
        // Declare inputs/outputs for dependency tracking.
        @Input(doc="unfixed bam") var unfixed: File = _
        def inputBams = List(unfixed)
        override def commandLine = super.commandLine + "%s%s%s IGNORE=INVALID_CIGAR IGNORE=MATE_NOT_FOUND".format(
          optional(" VALIDATION_STRINGENCY=", "SILENT"), repeat(" INPUT=", inputBams), " TMP_DIR=" + jobTempDir)
      }

      //fixMates.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
      validate.memoryLimit = 2
      validate.jarFile = qscript.picardValidateJar
      validate.unfixed = cleanedBam
      add(validate)

      val toQueryName = new PicardBamFunction {
        // Declare inputs/outputs for dependency tracking.
        @Input(doc="coordiante bam") var cobam: File = _
        @Output(doc="query bam") var qnbam: File = _
        def inputBams = List(cobam)
        def outputBam = qnbam
      }

      //fixMates.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
      toQueryName.memoryLimit = 4
      toQueryName.jarFile = qscript.picardSortSamJar
      toQueryName.cobam = cleanedBam
      toQueryName.qnbam = new File(cleanedBam.getAbsolutePath + ".qn.bam")
      add(toQueryName)

      Console.println("loop done")
    }
  }
}
