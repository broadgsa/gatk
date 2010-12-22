import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType
import scala.io.Source._

class PhaseSamples extends QScript {
  qscript => 
  
  @Input(doc="bam input, as .bam or as a list of files", shortName="I", required=true)
  var bams: File = _

  @Input(doc="master VCF calls", shortName="C", required=true)
  var masterCalls: File = _

  @Input(doc="file containing samples to be phased, as many as you like per line", shortName="samples", required=true)
  var samplesFile: File = _
  
  @Argument(doc="gatk jar file", shortName="J", required=true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc="ref", required=true)
  var referenceFile: File = _

  @Argument(fullName = "prefix", doc="Prefix argument", required=false)
  var prefix: String = ""

  @Argument(shortName = "L", doc="Intervals", required=false)
  var intervals: String = _

  @Input(doc="level of parallelism for BAM phaser.   By default is set to 0 [no scattering].", shortName="scatter", required=false)
  var scatterCount = 0

  @Input(doc="Samples to phase together.   By default is set to 1 [one job per sample].", shortName="samplesPerJob", required=false)
  var samplesPerJob = 1

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = Some(3)
    this.logging_level = "INFO"
  }

// A target has a list of samples and bam files to use for phasing
class Target(val name: String, val samples: List[String], val bams: List[File]) {
    def phasedVCFFile = new File(name + ".phased.vcf")
    
    override def toString(): String = String.format("[Target %s with samples %s against bams %s]", name, samples, bams)
}

def script = {
    if ( qscript.scatterCount > 0 ) throw new RuntimeException("scatter/gather currently not implemented")

    val samples: List[String] = getSamples(samplesFile)
    Console.out.printf("Samples are %s%n", samples)
    
    val targets: List[Target] = bamsToTargets(samples, bams)

    for (target <- targets) {
        Console.out.printf("Target is %s%n", target)
        add(new PhaseSamples(target))
    }

    add(new CombineVariants(targets.map(_.phasedVCFFile)))
}

def bamsToTargets(samples: List[String], bamsIn: File): List[Target] = {
    val bams: List[File] = parseBamsInput(bamsIn)
    
    def buildTargets(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
        case (Nil, y) => 
            return Nil 
        case (subsamples, remaining) => 
            return new Target("group" + count, subsamples, findBamsForSamples(subsamples, bams)) :: 
                    buildTargets(remaining, count + 1)
    }
    
    return buildTargets(samples, 0)
}    

// determines the bam files to use for these samples.  If there's only a single bam, assumes this file
// is a merge of all sample bams.  If there are multiple bams, for each sample, finds all bam files 
// containing sample in their path (should really check header?).  The phasing bams are the merge of 
// these files
def findBamsForSamples( samples: List[String], bams: List[File] ): List[File] = bams match {
    case x :: Nil => return bams
    case _ => 
        // todo -- implement me
        throw new RuntimeException("bam input lists not enabled")
}

// todo -- allow us to get the master list of samples from the VCF for convenience
def getSamples(samplesFile: File): List[String] = {
        return (for { line <- fromFile(samplesFile).getLines
                      part <- line.split("\\s", 0).iterator }
                        yield part ).toList
}

def parseBamsInput(bamsIn: File): List[File] = FilenameUtils.getExtension(bamsIn.getPath) match {
    case "bam" => return List(bamsIn)
    case "list" => return (for( line <- fromFile(bamsIn).getLines) yield new File(line)).toList
    case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam and .list")
}

class PhaseSamples(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.ReadBackedPhasing with CommandLineGATKArgs {
    this.rodBind :+= RodBind("variant", "VCF", qscript.masterCalls)
    this.out = t.phasedVCFFile
    this.input_file = t.bams
    this.sampleToPhase = t.samples
}

class CombineVariants(vcfsToCombine: List[File]) extends org.broadinstitute.sting.queue.extensions.gatk.CombineVariants with CommandLineGATKArgs {
    for ( vcf <- vcfsToCombine ) {
        this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    }

    // add the master call
    this.rodBind :+= RodBind("master", "VCF", masterCalls)
    this.variantMergeOptions = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.MASTER)

    this.out = new File("phased.vcf")
}

}
