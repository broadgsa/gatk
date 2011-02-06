package oneoffs.fromer

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

import org.broadinstitute.sting.gatk.DownsampleType


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

  @Output(doc = "DoC file to output", shortName = "o", required = true)
  var outputDoC: File = _

  @Input(doc = "Maximum depth (before GATK down-sampling kicks in...)", shortName = "MAX_DEPTH", required = false)
  var MAX_DEPTH = 20000

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    //this.memoryLimit = Some(3)
    this.logging_level = "INFO"
  }

  def script = {
    add(new DepthOfCoverage(bams, outputDoC))
  }

  class DepthOfCoverage(bam: File, docOutPath: File) extends org.broadinstitute.sting.queue.extensions.gatk.DepthOfCoverage with CommandLineGATKArgs {
    this.omitIntervalStatistics = false
    this.omitDepthOutputAtEachBase = true
    this.omitLocusTable = true

    this.minBaseQuality = Some(0)
    this.minMappingQuality = Some(0)

    this.out = docOutPath
    this.input_file :+= bam

    this.dcov = Some(MAX_DEPTH)
    this.downsampling_type = Some(DownsampleType.BY_SAMPLE)

    this.scatterCount = scatterCountInput

    override def dotString = "DOC: " + bam.getName
  }

}