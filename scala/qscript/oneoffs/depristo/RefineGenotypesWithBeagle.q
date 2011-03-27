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
  @Argument(doc="Memory in GB for beagle",required=false,shortName="BM") var BEAGLE_MEM_IN_GB: Int = 12
  @Argument(doc="X",required=false,shortName="cc") var CALIBRATION_CURVE: File = new File("vqsr.calibration.curve")

  @Argument(doc="X",required=false,shortName="test") var TEST: Boolean = false
  @Argument(doc="If provided, we'll skip over creating the reference panels, even if apparently required",required=false,shortName="assumeReferencePanelsExist") var assumeReferencePanelsExist: Boolean = false
  @Argument(doc="Tmpdir",required=false,shortName="tmpdir") var TMPDIR: File = new File("./")

  // assessing imputation performance
  @Argument(doc="VCF sites and alleles for genotyping assessment",required=false,shortName="assessmentSites") var assessmentSites: File = null
  @Argument(doc="BAM File for genotyping",required=false, shortName="bam") var bam: File = null
  @Argument(doc="Percent of sites that should be left out of BAM VCF to assess imputation", required=false, shortName="flo")
  var fractionsLeftOut: List[Double] = List(0.1, 0.2, 0.5, 0.9)
  // todo -- this might be best to think about in a different unit -- marker density per bp


  val MISSING_KEY = "?"
  val HM3_VCF: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf")
  val OMNI_VCF: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf")
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"

  trait GATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.jarFile = qscript.gatkJarFile
    this.memoryLimit = 2
  }

  // --------------------------------------------------------------------------------
  //
  // GENOTYPING SPECIFIC SITES IN A BAM FILE
  //
  // --------------------------------------------------------------------------------

  class GenotypeBAMAtSites(@Input bam: File, @Input sitesVCF: File, @Output genotypesVCF: File) extends UnifiedGenotyper with GATKArgs {
    this.input_file = List(bam)
    this.o = genotypesVCF
    this.stand_call_conf = 0.0
    this.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
    this.gt_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    this.rodBind :+= new RodBind("alleles","VCF",sitesVCF)

    // we only want chromosome counts annotations
    this.BTI = "alleles"
    this.G = List("none")
    this.A :+= "ChromosomeCounts"
    this.nsl = true

    // make sure we have the right intervals
    if ( interval != null ) {
      this.intervalsString = List(interval)
      this.BTIMR = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
    }
  }

  // --------------------------------------------------------------------------------
  //
  // BEAGLE COMMANDS
  //
  // --------------------------------------------------------------------------------

  class BeagleCommand(outputBase: String) extends CommandLineFunction {
    this.memoryLimit = BEAGLE_MEM_IN_GB

    // Note: These get set
    @Output val beaglePhasedFile: File = new File(outputBase +".phased.gz")
    @Output val beagleLikelihoods: File = new File(outputBase +".gprobs.gz")
    @Output val beagleRSquared: File = new File(outputBase +".r2")

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s out=ignore.me omitprefix=true".format(TMPDIR, BEAGLE_MEM_IN_GB, beagleJar)
  }

  class RefineGenotypesWithBeagle(@Input beagleInput: File, moreBeagleArgs: String = "")
    extends BeagleCommand(beagleInput.getName) {
    def myArgs = " like=%s %s".format(beagleInput.getAbsolutePath, moreBeagleArgs)
    override def commandLine = super.commandLine + myArgs
  }

  class ImputeMissingGenotypesWithReferencePanel(@Input evalBeagle: File,
                                                 @Input phasedBeagleFile: File,
                                                 @Input markers: File,
                                                 moreBeagleArgs: String = "")
    extends BeagleCommand(evalBeagle.getName) {
    def myArgs = " unphased=%s phased=%s markers=%s %s".format(evalBeagle.getAbsolutePath,
      phasedBeagleFile.getAbsolutePath, markers.getAbsolutePath, moreBeagleArgs)
    override def commandLine = super.commandLine + myArgs
  }

  class GunzipFile(@Input val in: File, @Output val out: File) extends CommandLineFunction {
    def commandLine = "gunzip -c %s > %s".format(in.getAbsolutePath, out.getAbsolutePath)
  }

  // --------------------------------------------------------------------------------
  //
  // CREATING AND EVALUATING REFERENCE PANELS
  //
  // --------------------------------------------------------------------------------

  class ReferencePanelBuilder(inputVCF: File, outputVCF: File, useCalibrationCurve: Boolean, moreBeagleArgs: String = "") {
    val beagleInput = new ProduceBeagleInput with GATKArgs
    if ( interval != null ) beagleInput.intervalsString = List(interval)
    beagleInput.variantVCF = inputVCF
    beagleInput.out = swapExt(outputVCF,".vcf",".beagle")
    if ( useCalibrationCurve ) beagleInput.cc = CALIBRATION_CURVE
    beagleInput.markers = swapExt(outputVCF, ".vcf", ".markers.txt")

    val refine = new RefineGenotypesWithBeagle(beagleInput.out, moreBeagleArgs)

    val unzipPhased = new GunzipFile(refine.beaglePhasedFile,swapExt(refine.beaglePhasedFile,".gz",".bgl"))
    val unzipProbs = new GunzipFile(refine.beagleLikelihoods,swapExt(refine.beagleLikelihoods,".gz",".bgl"))

    val vcfConvert = new BeagleOutputToVCF with GATKArgs
    vcfConvert.variantVCF = inputVCF
    vcfConvert.rodBind :+= new RodBind("beagleR2","BEAGLE",refine.beagleRSquared)
    vcfConvert.rodBind :+= new RodBind("beaglePhased","BEAGLE",unzipPhased.out)
    vcfConvert.rodBind :+= new RodBind("beagleProbs","BEAGLE",unzipProbs.out)
    vcfConvert.out = outputVCF

    def getMarkers = beagleInput.markers
    def getPanelPhasedHaplotypes = refine.beaglePhasedFile
    def getPanelVCF = vcfConvert.out

    def enqueueCommands() = {
      for ( cmd: CommandLineFunction <- List(beagleInput, refine, unzipPhased, unzipProbs, vcfConvert) )
        add(cmd)
    }
  }

  class EvaluateReferencePanel(@Input evalVCF: File,
                               @Output outputVCF: File,
                               panel: ReferencePanelBuilder,
                               percentLeftOut: Double,
                               moreBeagleArgs: String = "") {
    val evalBeagle = new VariantsToBeagleUnphased with GATKArgs
    if ( interval != null ) evalBeagle.intervalsString = List(interval)
    evalBeagle.variantVCF = evalVCF
    evalBeagle.out = swapExt(outputVCF,".vcf",".unphased.beagle")
    evalBeagle.bs = percentLeftOut
    evalBeagle.bsvcf = swapExt(outputVCF,".vcf",".missing.vcf")
    evalBeagle.missing = MISSING_KEY
    //evalBeagle.isIntermediate = true

    val refine = new ImputeMissingGenotypesWithReferencePanel(evalBeagle.out, panel.getPanelPhasedHaplotypes, panel.getMarkers, moreBeagleArgs)

    val unzipPhased = new GunzipFile(refine.beaglePhasedFile,swapExt(refine.beaglePhasedFile,".gz",".bgl"))
    val unzipProbs = new GunzipFile(refine.beagleLikelihoods,swapExt(refine.beagleLikelihoods,".gz",".bgl"))
    //unzipPhased.isIntermediate = true
    //unzipProbs.isIntermediate = true

    val vcfConvert = new BeagleOutputToVCF with GATKArgs
    vcfConvert.variantVCF = evalVCF
    vcfConvert.rodBind :+= new RodBind("beagleR2","BEAGLE",refine.beagleRSquared)
    vcfConvert.rodBind :+= new RodBind("beaglePhased","BEAGLE",unzipPhased.out)
    vcfConvert.rodBind :+= new RodBind("beagleProbs","BEAGLE",unzipProbs.out)
    vcfConvert.out = outputVCF
    vcfConvert.keep_monomorphic = true

    def getBootstrap: File = evalBeagle.bsvcf

    def enqueueCommands() = {
      for ( cmd: CommandLineFunction <- List(evalBeagle, refine, unzipPhased, unzipProbs, vcfConvert) )
        add(cmd)
    }
  }

  class EvalPanelAtChipSites(@Input eval: File) extends VariantEval with GATKArgs {
    this.noST = true
    this.evalModule :+= "GenotypeConcordance"
    this.o = swapExt(eval, ".vcf", ".vcf.eval")
    this.rodBind :+= RodBind("eval", "VCF", eval)
    this.rodBind :+= RodBind("comp_hm3", "VCF", HM3_VCF)
    this.rodBind :+= RodBind("comp_omni", "VCF", OMNI_VCF)
    if ( EvalInterval != null ) this.intervalsString = List(EvalInterval)
  }

  class EvalPanelAtBAMCalledSites(@Input imputedVCF: File, @Input bamGenotypes: File, @Input bootstrap: File) extends VariantEval with GATKArgs {
    this.evalModule :+= "GenotypeConcordance"
    this.o = swapExt(imputedVCF, ".vcf", ".vcf.eval")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37)
    this.rodBind :+= RodBind("eval", "VCF", imputedVCF)
    this.rodBind :+= RodBind("comp_bam_genotypes", "VCF", bamGenotypes)
    this.rodBind :+= RodBind("comp_bootstrap", "VCF", bootstrap)
    if ( EvalInterval != null ) this.intervalsString = List(EvalInterval)
  }

  def script = {
    var bamGenotypes: File = null

    if ( bam != null ) {
      bamGenotypes = new File(swapExt(bam, ".bam", "_genotyped_at." + assessmentSites.getName).getName)
      add(new GenotypeBAMAtSites(bam, assessmentSites, bamGenotypes))
    }

    for ( vcf <- vcfsToBeagle ) {
      if ( ! TEST ) add(new EvalPanelAtChipSites(vcf))
      for ( useCalibrationCurve <- List(true, false) ) {
        for ( niter <- List(10, 20, 50) ) {
          if ( ! TEST || (niter == 10 && ! useCalibrationCurve)) {
            val refineFilenamePart = "niter_%d_cc_%b".format(niter, useCalibrationCurve)
            val refine_out = swapExt(vcf, ".vcf", ".refined.%s.vcf".format(refineFilenamePart))
            val refPanel = new ReferencePanelBuilder(vcf, refine_out, useCalibrationCurve, "niterations=%d".format(niter))
            if ( ! assumeReferencePanelsExist ) refPanel.enqueueCommands()

            // start up VE
            add(new EvalPanelAtChipSites(refine_out))

            if ( bamGenotypes != null ) {
              for ( fractionLeftOut <- if ( TEST ) List(0.1) else fractionsLeftOut ) {
                val bamGenotypesImputed = swapExt(bamGenotypes, ".vcf", "_flo_%.2f.imputed_with_%s".format(fractionLeftOut, refine_out))
                val args = "niterations=%d missing=\"%s\"".format(niter, MISSING_KEY)
                val panelEval = new EvaluateReferencePanel(bamGenotypes, bamGenotypesImputed, refPanel, fractionLeftOut, args)
                panelEval.enqueueCommands()
                add(new EvalPanelAtBAMCalledSites(bamGenotypesImputed, bamGenotypes, panelEval.getBootstrap))
              }
            }
          }
        }
      }
    }
  }
}