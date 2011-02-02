import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils
import tools.nsc.io.Process;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType
import scala.io.Source._

class PhaseSamples extends QScript {
  qscript =>

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "master VCF calls", shortName = "C", required = true)
  var masterCalls: File = _

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(fullName = "prefix", doc = "Prefix argument", required = false)
  var prefix: String = ""

  @Argument(shortName = "L", doc = "Intervals", required = false)
  var intervals: String = _

  @Input(doc = "level of parallelism for BAM phaser.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCount = 0

  @Input(doc = "Samples to phase together.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Phased file to output", shortName = "o", required = true)
  var outputPhased: File = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = Some(3)
    this.logging_level = "INFO"
  }

  // A target has a list of samples and bam files to use for phasing
  class Target(val name: String, val samples: List[String], val bams: List[File]) {
    def phasedVCFFile = new File(name + "." + outputPhased)

    override def toString(): String = String.format("[Target %s with samples %s against bams %s]", name, samples, bams)
  }

  def script = {
    if (qscript.scatterCount > 0) throw new RuntimeException("scatter/gather currently not implemented")

    val samples: List[String] = getSamples()
    Console.out.printf("Samples are %s%n", samples)

    val targets: List[Target] = bamsToTargets(samples, bams)

    for (target <- targets) {
      Console.out.printf("Target is %s%n", target)
      add(new PhaseSamples(target))
    }

    add(new CombineVariants(targets.map(_.phasedVCFFile)))

    add(new PhasingByACeval())
  }

  def bamsToTargets(samples: List[String], bamsIn: File): List[Target] = {
    val bams: List[File] = parseBamsInput(bamsIn)
    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = getMapOfBamsForSample(bams)

    def buildTargets(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
      case (Nil, y) =>
        return Nil
      case (subsamples, remaining) =>
        return new Target("group" + count, subsamples, findBamsForSamples(subsamples, sampleToBams)) ::
                buildTargets(remaining, count + 1)
    }

    return buildTargets(samples, 0)
  }

  def getMapOfBamsForSample(bams: List[File]): scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = bams match {
    case Nil => return scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]

    case x :: y =>
      val m: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = getMapOfBamsForSample(y)

      val getBamSamplesCommand: String = "samtools view -H " + x.getPath() + " | grep '^@RG' | awk '{for (i = 1; i <= NF; i++) if (substr($i,1,3) == \"SM:\") print substr($i,4)}' | sort | uniq"
      //println("getBamSamplesCommand: " + getBamSamplesCommand)
      val bamSamples: List[String] = Process(getBamSamplesCommand).iterator toList

      for (s <- bamSamples) {
        if (!m.contains(s))
          m += s -> scala.collection.mutable.Set.empty[File]

        m(s) = m(s) + x
      }

      return m
  }

  def findBamsForSamples(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[File] = {
    val l: List[File] = Nil
    l ++ findBamsForSamplesHelper(samples, sampleToBams)
  }

  def findBamsForSamplesHelper(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): scala.collection.mutable.Set[File] = samples match {
    case Nil => scala.collection.mutable.Set.empty[File]

    case x :: y =>
      var bamsForSampleX: scala.collection.mutable.Set[File] = scala.collection.mutable.Set.empty[File]
      if (sampleToBams.contains(x))
        bamsForSampleX = sampleToBams(x)
      return bamsForSampleX ++ findBamsForSamplesHelper(y, sampleToBams)
  }

  def getSamples(): List[String] = {
    val getSamplesCommand: String = "cat " + masterCalls.getPath() + " | grep '^#CHROM' | head -1 | awk '{for (i = 10; i <= NF; i++) print $i}'"
    //println("getSamplesCommand: " + getSamplesCommand)
    return Process(getSamplesCommand).iterator toList
  }

  def parseBamsInput(bamsIn: File): List[File] = FilenameUtils.getExtension(bamsIn.getPath) match {
    case "bam" => return List(bamsIn)
    case "list" => return (for (line <- fromFile(bamsIn).getLines) yield new File(line)).toList
    case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam and .list")
  }

  class PhaseSamples(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.ReadBackedPhasing with CommandLineGATKArgs {
    this.rodBind :+= RodBind("variant", "VCF", qscript.masterCalls)
    this.out = t.phasedVCFFile
    this.input_file = t.bams
    this.sampleToPhase = t.samples
  }

  class CombineVariants(vcfsToCombine: List[File]) extends org.broadinstitute.sting.queue.extensions.gatk.CombineVariants with CommandLineGATKArgs {
    for (vcf <- vcfsToCombine) {
      this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    }

    // add the master call:
    this.rodBind :+= RodBind("master", "VCF", masterCalls)
    this.variantMergeOptions = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.MASTER)

    this.out = outputPhased
  }

  class PhasingByACeval() extends org.broadinstitute.sting.queue.extensions.gatk.PhasingEval with CommandLineGATKArgs {
    this.analysis = org.broadinstitute.sting.oneoffprojects.walkers.phasing.PhasingEval.Analysis.PHASING_BY_AC

    this.rodBind :+= RodBind(outputPhased.getName, "VCF", outputPhased)

    this.out = new File("phasing_by_ac." + outputPhased + ".txt")
  }

}