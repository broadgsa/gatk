package oneoffs.fromer.CNV

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import java.io.PrintWriter
import org.apache.commons.io.IOUtils

class ReadDepthCNVanalysis extends QScript {
  qscript =>

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(shortName = "L", doc = "Intervals", required = false)
  var intervals: String = _

  @Input(doc = "level of parallelism for BAM DoC.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCountInput = 0

  @Input(doc = "Samples to phase together.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "DoC file to output", shortName = "o", required = true)
  var outputDoC: File = _

  @Input(doc = "Maximum depth (before GATK down-sampling kicks in...)", shortName = "MAX_DEPTH", required = false)
  var MAX_DEPTH = 20000

  @Input(doc = "Number of read-depth bins", shortName = "NUM_BINS", required = false)
  var NUM_BINS = 200

  @Input(doc = "Starting value of read-depth bins", shortName = "START_BIN", required = false)
  var START_BIN = 1

  @Input(doc = "Minimum read mapping quality", shortName = "MMQ", required = false)
  var minMappingQuality = 0

  val DOC_OUTPUT_SUFFIX: String = ".sample_interval_summary"

  val DOC_MEAN_COVERAGE_OUTPUT: String = ".sample_interval.averageCoverage.txt"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    //this.memoryLimit = 3
    this.logging_level = "INFO"
  }

  // A target has a list of samples and bam files to use for DoC
  class Target(val name: String, val samples: List[String], val bams: List[File]) {
    var prefix: String = outputDoC.getParent()
    if (prefix == null)
      prefix = ""
    else
      prefix = prefix + "/"

    def DoC_output = new File(prefix + name + "." + outputDoC.getName())

    override def toString(): String = String.format("[Target %s [%s] with samples %s against bams %s]", name, DoC_output, samples, bams)
  }

  def script = {
    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = VCF_BAM_utilities.getMapOfBAMsForSample(VCF_BAM_utilities.parseBAMsInput(bams))
    val samples: List[String] = sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val targets: List[Target] = buildTargets(samples, sampleToBams)

    for (target <- targets) {
      Console.out.printf("Target is %s%n", target)
      add(new DoC(target))
    }

    add(new combineDoC(targets.map(u => new File(u.DoC_output.getPath() + DOC_OUTPUT_SUFFIX))))
  }

  def buildTargets(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[Target] = {

    def buildTargetsHelper(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
      case (Nil, y) =>
        return Nil
      case (subsamples, remaining) =>
        return new Target("group" + count, subsamples, VCF_BAM_utilities.findBAMsForSamples(subsamples, sampleToBams)) ::
          buildTargetsHelper(remaining, count + 1)
    }

    return buildTargetsHelper(samples, 0)
  }

  class DoC(t: Target) extends CommandLineGATKArgs with ScatterGatherableFunction {
    this.analysis_type = "DepthOfCoverage"

    this.input_file = t.bams

    this.downsample_to_coverage = MAX_DEPTH
    this.downsampling_type = DownsampleType.BY_SAMPLE

    this.scatterCount = scatterCountInput
    this.scatterClass = classOf[IntervalScatterFunction]

    @Output
    @Gather(classOf[org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction])
    var intervalSampleOut: File = new File(t.DoC_output.getPath() + DOC_OUTPUT_SUFFIX)

    val outFile = new File(intervalSampleOut.getParentFile(), t.DoC_output.getName())

    override def commandLine = super.commandLine +
      " --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality " + minMappingQuality +
      " --start " + START_BIN + " --stop " + MAX_DEPTH + " --nBins " + NUM_BINS +
      " -o " + outFile

    override def dotString = "DOC: " + t.DoC_output

    this.jobOutputFile = outFile + ".out"
  }

  class combineDoC(DoCsToCombine: List[File]) extends CommandLineFunction {
    override def description = "Combines DoC outputs for multiple samples (at same loci)"

    @Input(doc = "")
    var inputDoCfiles: List[File] = DoCsToCombine

    @Output
    val outputDoCaverageCoverage: File = new File(outputDoC.getPath + DOC_MEAN_COVERAGE_OUTPUT)

    var command: String = "~fromer/CNV/wave1+2/scripts/mergeDoC.pl -gatk " + qscript.gatkJarFile.getPath.replaceFirst("dist/GenomeAnalysisTK.jar", "") + " -ref " + qscript.referenceFile + " -out " + outputDoCaverageCoverage
    for (input <- inputDoCfiles) {
      command += " " + input
    }
    def commandLine = command

    // Since loading ALL of the output into the perl script can take significant memory:
    this.memoryLimit = 9
  }
}