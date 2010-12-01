import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK
import org.broadinstitute.sting.queue.pipeline.{BamProcessing,VariantCalling}
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils

class fullCallingPipelineV2 extends QScript {
  qscript =>

  @Argument(doc="Number of cleaning jobs", shortName="cleaningJobs", required=false)
  var cleaningJobs: Int = 1

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="path to trigger track (for UnifiedGenotyper)", shortName="trigger", required=false)
  var trigger: File = _

  @Input(doc="path to refseqTable (for GenomicAnnotator)", shortName="refseqTable")
  var refseqTable: File = _

  @Input(doc="path to Picard FixMateInformation.jar.  See http://picard.sourceforge.net/ .", required=false)
  var picardFixMatesJar: File = new java.io.File("/seq/software/picard/current/bin/FixMateInformation.jar")

  @Input(doc="path to GATK jar", shortName="gatk")
  var gatkJar: File = _

  @Input(doc="target Ti/Tv ratio for recalibration", shortName="titv", required=true)
  var target_titv: Float = _

  @Input(doc="per-sample downsampling level",shortName="dcov",required=false)
  var downsampling_coverage = 300

  @Input(doc="level of parallelism for UnifiedGenotyper", shortName="snpScatter", required=false)
  var num_snp_scatter_jobs = 20

  @Input(doc="level of parallelism for IndelGenotyperV2", shortName="indelScatter", required=false)
  var num_indel_scatter_jobs = 5

  @Input(doc="Skip indel-cleaning for BAM files (for testing only)", shortName="skipCleaning", required=false)
  var skip_cleaning = false

  //@Input(doc="ADPR script")
  //var adprScript: File = _

  //@Input(doc="Sequencing maching name (for use by adpr)")
  //var machine: String = _

  //@Input(doc="Sequencing experiement type (for use by adpr)--Whole_Exome, Whole_Genome, or Hybrid_Selection")
  //var protocol: String = _

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals :+= qscript.pipeline.getProject.getIntervalList
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.memoryLimit = Some(4)
  }


  // ------------ SETUP THE PIPELINE ----------- //


  def script = {
    pipeline = YamlUtils.load(classOf[Pipeline], qscript.yamlFile)
    var callingLib: VariantCalling = new VariantCalling(qscript.yamlFile,qscript.gatkJar)
    var cleaningLib: BamProcessing = new BamProcessing(qscript.yamlFile,qscript.gatkJar,qscript.picardFixMatesJar)

    val projectBase: String = qscript.pipeline.getProject.getName
    val cleanedBase: String = projectBase + ".cleaned"
    val uncleanedBase: String = projectBase + ".uncleaned"

    // there are commands that use all the bam files
    val recalibratedSamples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains("recalibrated"))

    var bamsToClean: List[(File,File)] = Nil
    var recalBams: List[File] = Nil
    var cleanedBams: List[File] = Nil

    for ( sample <- recalibratedSamples ) {
      val bam = sample.getBamFiles.get("recalibrated")
      recalBams :+= bam
      if (!sample.getBamFiles.contains("cleaned")) {
        sample.getBamFiles.put("cleaned", swapExt(bam,"bam","cleaned.bam"))
        bamsToClean :+= (bam,sample.getBamFiles.get("cleaned"))
      }

      cleanedBams :+= sample.getBamFiles.get("cleaned")
    }

    if ( !qscript.skip_cleaning ) {
      addAll(cleaningLib.StandardIndelRealign(bamsToClean,qscript.cleaningJobs))
    }

    if (!qscript.skip_cleaning) {
      endToEnd(cleanedBase, cleanedBams, callingLib)
    } else {
      endToEnd(uncleanedBase, recalBams, callingLib)
    }
  }


  def endToEnd(base: String, bamFiles: List[File], lib: VariantCalling) = {
    var recal_vcf = new File(base+"_snps.recal.annotated.tranched.vcf")
    var handfilt_vcf = new File(base+"_snps.handfiltered.annotated.vcf")
    var indel_vcf = new File(base+"_indel_calls.vcf")

    addAll(lib.StandardCallingPipeline(bamFiles,indel_vcf,recal_vcf,handfilt_vcf,qscript.target_titv,qscript.refseqTable))
  }

  def addAll(clfs: List[CommandLineFunction]) = {
    for ( c <- clfs ) {
      add(c)
    }
  }
}
