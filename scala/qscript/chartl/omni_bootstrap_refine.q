import java.io.{FileReader, File, BufferedReader}
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

class omni_bootstrap_refine extends QScript {
 => script

  val ref = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val gatkJar = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")

  trait FooTrait extends CommandLineGATK {
    this.reference_sequence = script.ref
    this.jarFile = script.gatkJar
    this.intervalsString :+= "20"
  }

  class SampleOverlap extends CommandLineFunction {
    @Input(doc="omni vcf") var omni: File = _
    @Input(doc="other vcf") var other: File = _
    @Output(doc="output sn file") var out : File = _

    def commandLine = "%s ; %s ; %s ; %s".format(this.cmd1, this.cmd2, this.cmd3, this.cmd4)

    def cmd1 : String = {
      return "head -n 500 %s | grep #CHR | cut -f10- | tr '\t' '\n' | sort > temp1".format(omni.getAbsolutePath)
    }

    def cmd2 : String = {
      return "head -n 500 %s | grep #CHR | cut -f10- | tr '\t' '\n'| sort > temp2".format(other.getAbsolutePath) 
    }

    def cmd3 : String = {
      return "cat temp1 temp2 | sort | uniq -c | awk '{if($1 == 2) print $2}' > %s".format(out.getAbsolutePath)
    }

    def cmd4 : String = {
      return "rm temp1 temp2"
    }
  }

  class BeagleRefinement extends CommandLineFunction {
    @Input(doc="The beagle input file") var beagleInput: File = _
    var beagleOutputBase: String = _
    var beagleMemoryGigs: Int = 4

    /**
     * Note: These get set
     */
    @Output(doc="The beagle phased file") var beaglePhasedFile: File = _
    @Output(doc="The beagle likelihood file") var beagleLikelihoods: File = _
    @Output(doc="The beagle r2 file") var beagleRSquared: File = _
    var beagleOutputDir: String = _

    def freezeOutputs = {
      if ( beagleInput.getParent == null ) {
        beagleOutputDir = ""
      } else {
        beagleOutputDir = beagleInput.getParent
      }
      beaglePhasedFile = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".phased.gz")
      beagleLikelihoods = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".gprobs.gz")
      beagleRSquared = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".r2")
    }

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s like=%s out=%s".format(beagleInput.getParent,beagleMemoryGigs,beagleJar,beagleInput.getAbsolutePath,beagleOutputBase)
  }

  def runme(pop: File, omni: File, base: String) : List[CommandLineFunction] = {
    var clf : List[CommandLineFunction] = Nil
    val s_overlap = new File(base+"_sample_overlap.txt")
    var calcOverlap = new SampleOverlap
    calcOverlap.omni = omni
    calcOverlap.other = pop
    calcOverlap.out = s_overlap

    clf += calcOverlap

    val subset_omni = swapExt(omni,".vcf","_%s_subset.vcf".format(base))
    val subset_pop = swapExt(pop,".vcf","_%s_subset.vcf".format(base))

    var omSubset = new SelectVariants with FooTrait
    omSubset.variantVCF = omni
    omSubset.sample = s_overlap.getAbsolutePath
    omSubset.out = subset_omni

    clf += omSubset

    var popSubset = new SelectVariants with FooTrait
    popSubset.variantVCF = pop
    popSubset.sample = s_overlap.getAbsolutePath
    popSubset.out = subset_pop

    clf += popSubset

    var bootRefine = new ProduceBeagleInput with FooTrait
    bootRefine.variantVCF = popSubset
    bootRefine.rodBind :+= new RodBind("validation","VCF",subset_omni)
    bootRefine.bootstrap_vcf = swapExt(subset_omni,".vcf",".boot.vcf")
    bootRefine.bootstrap = Some(2)
    bootRefine.validation_genotype_ptrue = Some(0.98);
    bootRefine.out = new File(base+".beagle_input.beagle")

    clf += bootRefine

    var runBeagle = new BeagleRefinement
    runBeagle.beagleInput = bootRefine.out
    runBeagle.beagleOutputBase = base+"_beagle"
    runBeagle.freezeOutputs()

    clf += runBeagle

    var putBack = new BeagleOutputToVCF with FooTrait
    putBack.
  }
}