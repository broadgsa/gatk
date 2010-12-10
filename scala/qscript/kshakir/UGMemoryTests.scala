import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import collection.JavaConversions._

class UGMemoryTests extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="The path to the GenomeAnalysisTK.jar file.", shortName="gatk")
  var gatkJar: File = null

  @Input(doc="per-sample downsampling level",shortName="dcov",required=false)
  var downsampling_coverage = 300

  def script = {
    val pipeline = YamlUtils.load(classOf[Pipeline], qscript.yamlFile)
    val memoryLimits = List(1,2,4,6,8,10,12,16)
    val recalibratedSamples = pipeline.getSamples.map(_.getBamFiles.get("recalibrated")).toList
    val squid1 = "C315"
    val squid2 = "C338"
    
    val numBamsList = List(10, 20, 50, 70, 100, 120, 150)
    val squid1Bams = recalibratedSamples.filter(_.getAbsolutePath.contains(squid1))
    val squid2Bams = recalibratedSamples.filter(_.getAbsolutePath.contains(squid2))

    for (memoryLimit <- memoryLimits) {
      for (numBams <- numBamsList) {
        val dir = "%03d_bams_%02dg".format(numBams, memoryLimit)

        val snps = new UnifiedGenotyper
        snps.jobOutputFile = new File(dir, "UnifiedGenotyper.out")
        snps.out = new File(dir, "UnifiedGenotyper.vcf")
        snps.input_file = squid1Bams.take(numBams/2) ++ squid2Bams.take(numBams/2)
        snps.memoryLimit = Some(memoryLimit)

        snps.jarFile = qscript.gatkJar
        snps.reference_sequence = pipeline.getProject.getReferenceFile
        snps.intervals = List(pipeline.getProject.getIntervalList)
        snps.DBSNP = pipeline.getProject.getDbsnpFile
        snps.downsample_to_coverage = Some(qscript.downsampling_coverage)
        snps.annotation ++= List("AlleleBalance")
        snps.group :+= "Standard"

        add(snps)
      }
    }
  }
}
