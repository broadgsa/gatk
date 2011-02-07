package oneoffs.fromer

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.util.BAMutilities


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

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    //this.memoryLimit = Some(3)
    this.logging_level = "INFO"
  }

  // A target has a list of samples and bam files to use for DoC
  class Target(val name: String, val samples: List[String], val bams: List[File]) {
    def DoC_output = new File(name + "." + outputDoC)

    override def toString(): String = String.format("[Target %s with samples %s against bams %s]", name, samples, bams)
  }

  def script = {
    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = BAMutilities.getMapOfBamsForSample(BAMutilities.parseBamsInput(bams))
    val samples: List[String] = sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val targets: List[Target] = buildTargets(samples, sampleToBams)

    for (target <- targets) {
      Console.out.printf("Target is %s%n", target)
      add(new DoC(target))
    }
  }

  def buildTargets(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[Target] = {

    def buildTargetsHelper(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
      case (Nil, y) =>
        return Nil
      case (subsamples, remaining) =>
        return new Target("group" + count, subsamples, BAMutilities.findBamsForSamples(subsamples, sampleToBams)) ::
                buildTargetsHelper(remaining, count + 1)
    }

    return buildTargetsHelper(samples, 0)
  }

  class DoC(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.DepthOfCoverage with CommandLineGATKArgs {
    this.omitIntervalStatistics = false
    this.omitDepthOutputAtEachBase = true
    this.omitLocusTable = true

    this.minBaseQuality = Some(0)
    this.minMappingQuality = Some(0)

    this.out = t.DoC_output
    this.input_file = t.bams

    this.dcov = Some(MAX_DEPTH)
    this.downsampling_type = Some(DownsampleType.BY_SAMPLE)

    this.start = Some(START_BIN)
    this.stop = Some(MAX_DEPTH)
    this.nBins = Some(NUM_BINS)

    this.scatterCount = scatterCountInput

    override def dotString = "DOC: " + t.DoC_output
  }

}