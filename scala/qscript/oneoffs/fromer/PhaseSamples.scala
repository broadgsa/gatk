package oneoffs.fromer

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities

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
  var intervals: String = null

  @Input(doc = "level of parallelism for BAM phaser.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCount = 0

  @Input(doc = "Samples to phase together.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Phased file to output", shortName = "o", required = true)
  var outputPhased: File = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    if (qscript.intervals != null) {
      this.intervalsString = List(qscript.intervals)
    }
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = 3
    this.logging_level = "INFO"
  }

  // A target has a list of samples and bam files to use for phasing
  class Target(val name: String, val samples: List[String], val bams: List[File]) {
    var prefix: String = outputPhased.getParent()
    if (prefix == null)
      prefix = ""
    else
      prefix = prefix + "/"

    def phasedVCFFile = new File(prefix + name + "." + outputPhased.getName())

    override def toString(): String = String.format("[Target %s [%s] with samples %s against bams %s]", name, phasedVCFFile, samples, bams)
  }

  def script = {
    if (qscript.scatterCount > 0) throw new RuntimeException("scatter/gather currently not implemented")

    val samples: List[String] = VCF_BAM_utilities.getSamplesFromVCF(masterCalls)
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
    val bams: List[File] = VCF_BAM_utilities.parseBAMsInput(bamsIn)
    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = VCF_BAM_utilities.getMapOfBAMsForSample(bams)

    def buildTargets(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
      case (Nil, y) =>
        return Nil
      case (subsamples, remaining) =>
        return new Target("group" + count, subsamples, VCF_BAM_utilities.findBAMsForSamples(subsamples, sampleToBams)) ::
                buildTargets(remaining, count + 1)
    }

    return buildTargets(samples, 0)
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
    this.variantMergeOptions = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.MASTER

    this.out = outputPhased
  }

  class PhasingByACeval() extends org.broadinstitute.sting.queue.extensions.gatk.PhasingEval with CommandLineGATKArgs {
    this.analysis = org.broadinstitute.sting.oneoffprojects.walkers.phasing.PhasingEval.Analysis.PHASING_BY_AC

    this.rodBind :+= RodBind(outputPhased.getName, "VCF", outputPhased)

    this.out = new File("phasing_by_ac." + outputPhased + ".txt")
  }

}