package oneoffs.depristo

//import net.sf.picard.reference.FastaSequenceFile
//import org.broadinstitute.sting.datasources.pipeline.Pipeline
//import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
//import collection.JavaConversions._

class RefineGenotypesWithBeagle extends QScript {
  qscript =>

  @Argument(doc="VCF file to run beagle genotype refinement on",required=true,shortName="vcf") var vcfsToBeagle: List[File] = _
  @Argument(doc="Path to GATK jar",required=true,shortName="gatkjarfile") var gatkJarFile: File = _
  @Argument(doc="Path to BEAGLE jar",required=true,shortName="beagle") var beagleJar: File = _
  @Argument(doc="Reference file",required=true,shortName="R") var reference: File = _
  @Argument(doc="Beagle interval",required=false,shortName="L") var interval: String = null
  @Argument(doc="Evaluation interval",required=false,shortName="Le") var EvalInterval: String = null
  @Argument(doc="Memory in GB for beagle",required=false,shortName="BM") var BEAGLE_MEM_IN_GB: Int = 6
  @Argument(doc="X",required=false,shortName="cc") var CALIBRATION_CURVE: File = new File("vqsr.calibration.curve")

  @Argument(doc="X",required=false,shortName="test") var TEST: Boolean = false

  val HM3_VCF: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf")
  val OMNI_VCF: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf")

  trait GATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.jarFile = qscript.gatkJarFile
    this.memoryLimit = Some(2)
  }

  class GunzipFile(in: File, out:File ) extends CommandLineFunction {
    @Input(doc="file to gunzip") var inp = in
    @Output(doc="file to gunzip to") var outp = out

    def commandLine = "gunzip -c %s > %s".format(inp.getAbsolutePath, outp.getAbsolutePath)
  }

  class BeagleRefinement(moreBeagleArgs: String = "") extends CommandLineFunction {
    @Input(doc="The beagle input file") var beagleInput: File = _
    var beagleOutputBase: String = _
    var beagleMemoryGigs: Int = BEAGLE_MEM_IN_GB

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
      beaglePhasedFile = new File(beagleOutputDir + beagleOutputBase+"."+beagleInput.getName+".phased.gz")
      beagleLikelihoods = new File(beagleOutputDir + beagleOutputBase+"."+beagleInput.getName+".gprobs.gz")
      beagleRSquared = new File(beagleOutputDir + beagleOutputBase +"."+beagleInput.getName+".r2")
    }

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s like=%s %s out=%s".format(beagleInput.getParent,beagleMemoryGigs,beagleJar,beagleInput.getAbsolutePath,moreBeagleArgs,beagleOutputBase)
  }

  def RefineGenotypes(inputVCF: File, outputVCF: File, useCalibrationCurve: Boolean, moreBeagleArgs: String = "") = {
    var beagleInput = new ProduceBeagleInput with GATKArgs
    if ( interval != null ) beagleInput.intervalsString = List(interval)
    beagleInput.variantVCF = inputVCF
    beagleInput.out = swapExt(outputVCF,".vcf",".beagle")
    if ( useCalibrationCurve ) beagleInput.cc = CALIBRATION_CURVE

    var refine = new BeagleRefinement(moreBeagleArgs)
    refine.beagleInput = beagleInput.out
    refine.beagleOutputBase = outputVCF.getName + ".bout"
    refine.beagleMemoryGigs = BEAGLE_MEM_IN_GB
    refine.memoryLimit = Some(BEAGLE_MEM_IN_GB)
    refine.freezeOutputs

    var unzipPhased = new GunzipFile(refine.beaglePhasedFile,swapExt(refine.beaglePhasedFile,".gz",".bgl"))
    var unzipProbs = new GunzipFile(refine.beagleLikelihoods,swapExt(refine.beagleLikelihoods,".gz",".bgl"))
//    unzipPhased.isIntermediate = true
//    unzipProbs.isIntermediate = true

    var vcfConvert = new BeagleOutputToVCF with GATKArgs
    vcfConvert.variantVCF = inputVCF
    vcfConvert.rodBind :+= new RodBind("beagleR2","BEAGLE",refine.beagleRSquared)
    vcfConvert.rodBind :+= new RodBind("beaglePhased","BEAGLE",unzipPhased.outp)
    vcfConvert.rodBind :+= new RodBind("beagleProbs","BEAGLE",unzipProbs.outp)
    vcfConvert.out = outputVCF

    for ( cmd: CommandLineFunction <- List(beagleInput, refine, unzipPhased, unzipProbs, vcfConvert) )
      add(cmd)
  }

  class MyEval(@Input(doc="Evaluation file") eval: File) extends VariantEval with GATKArgs {
    this.noST = true
    this.evalModule :+= "GenotypeConcordance"
    this.o = swapExt(eval, ".vcf", ".vcf.eval")
    this.rodBind :+= RodBind("eval", "VCF", eval)
    this.rodBind :+= RodBind("comp_hm3", "VCF", HM3_VCF)
    this.rodBind :+= RodBind("comp_omni", "VCF", OMNI_VCF)
    if ( EvalInterval != null ) {
      Console.printf("EvalInterval " + EvalInterval)
      this.intervalsString = List(EvalInterval)
    }
  }

  def script = {
    for ( vcf <- vcfsToBeagle ) {
      if ( ! TEST ) add(new MyEval(vcf))
//      for ( niter <- List(10) ) {
      for ( useCalibrationCurve <- List(true, false) ) {
//        for ( includeFilteredRecords <- List(true, false) ) {
        for ( niter <- List(10, 25, 50) ) {
          if ( ! TEST || (niter == 10 && ! useCalibrationCurve)) {
            val refine_out = swapExt(vcf, ".vcf", ".refined.niter_%d_cc_%b.vcf".format(niter, useCalibrationCurve))
            RefineGenotypes(vcf, refine_out, useCalibrationCurve, "niterations=%d".format(niter))
            add(new MyEval(refine_out))
          }
        }
      }
    }
  }
}