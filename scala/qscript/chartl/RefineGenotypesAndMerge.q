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

class BeagleGenotypeRefinement extends QScript {
  qscript =>

  @Argument(doc="VCF file to run beagle genotype refinement on",required=true,shortName="vcf") var vcfsToBeagle: List[File] = _
  @Argument(doc="Output file to which to write final vcf",required=true,shortName="out") var outVCF: File = _
  @Argument(doc="Path to GATK jar",required=true,shortName="gatk") var gatkJar: File = _
  @Argument(doc="Path to BEAGLE jar",required=true,shortName="beagle") var beagleJar: File = _
  @Argument(doc="Reference file",required=true,shortName="ref") var reference: File = _

  trait GATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.jarFile = qscript.gatkJar
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

    def freezeOutputs = {
      val beagleOutputDir = beagleInput.getParent
      beaglePhasedFile = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".phased")
      beagleLikelihoods = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".gprobs")
      beagleRSquared = new File(beagleOutputDir+beagleOutputBase+"."+beagleInput.getName+".r2")
    }

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s like=%s out=%s".format(beagleInput.getParent,beagleMemoryGigs,beagleJar,beagleInput.getAbsolutePath,beagleOutputBase)
  }

  def RefineGenotypes(inputVCF: File, outputVCF: File, beagleBase: String ) : List[CommandLineFunction] = {
    var commands: List[CommandLineFunction] = Nil

    var beagleInput = new ProduceBeagleInput with GATKArgs
    beagleInput.variantVCF = inputVCF
    beagleInput.out = swapExt(inputVCF,".vcf",".beagle")

    var refine = new BeagleRefinement
    refine.beagleInput = beagleInput.out
    refine.beagleOutputBase = beagleBase
    refine.beagleMemoryGigs = 20
    refine.memoryLimit = Some(4)
    refine.freezeOutputs

    var vcfConvert = new BeagleOutputToVCF with GATKArgs
    vcfConvert.variantVCF = inputVCF
    vcfConvert.rodBind :+= new RodBind("beagleR2","BEAGLE",refine.beagleRSquared)
    vcfConvert.rodBind :+= new RodBind("beaglePhased","BEAGLE",refine.beaglePhasedFile)
    vcfConvert.rodBind :+= new RodBind("beagleProbs","BEAGLE",refine.beagleLikelihoods)
    vcfConvert.out = outputVCF

    commands :+= beagleInput
    commands :+= refine
    commands :+= vcfConvert

    return commands
  }

  def mergeVCFs(vcfs: List[File], outputVCF: File) : CombineVariants = {
    var cv = new CombineVariants with GATKArgs
    cv.out = outputVCF
    cv.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    cv.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    cv.priority = (vcfs.foldLeft[List[String]](Nil)( (bNames,vcf) => bNames ::: List[String](swapExt(vcf,".vcf","").getName))).mkString(",")
    cv.rodBind = vcfs.foldLeft[List[RodBind]](Nil)( (rods,vcf) => rods ::: List[RodBind](new RodBind(swapExt(vcf,".vcf","").getName,"VCF",vcf)))

    return cv
  }

  def script = {
    var vcfsToMerge: List[File] = Nil
    for ( i <- 0 until vcfsToBeagle.size ) {
      val base = vcfsToBeagle.get(i).getName+".bout"
      val refine_out = swapExt(vcfsToBeagle.get(i),".vcf",".refined.vcf")
      vcfsToMerge :+= refine_out
      for ( c <- RefineGenotypes(vcfsToBeagle.get(i),refine_out,base) ) {
        add(c)
      }
    }

    add(mergeVCFs(vcfsToMerge,outVCF))
  }
}