package oneoffs.depristo

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class AssessChipCoverageOfPanel extends QScript {
  qscript =>
  @Argument(doc="Path to GATK jar",required=false,shortName="gatkjarfile") var gatkJarFile: File = new File("dist/GenomeAnalysisTK.jar")
  @Argument(doc="Panel VCF",required=true,shortName="panelVCF") var panelVCF: File = _
  @Argument(doc="BAM File",required=true, shortName="bam") var bam: File = null
  @Argument(doc="Bundle path",required=false, shortName="bundle") var bundle: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")

  @Argument(shortName = "R", doc="ref", required=true)
  var referenceFile: File = _

  @Argument(shortName = "L", doc="intervals", required=false)
  val TARGET_INTERVAL: String = null;

  def HM3_VCF: File = new File(bundle + "/hapmap_3.3.b37.sites.vcf")
  def OMNI_VCF: File = new File(bundle + "/1000G_omni2.5.b37.sites.vcf")

  trait GATKArgs extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    if ( TARGET_INTERVAL != null )
      this.intervalsString = List(TARGET_INTERVAL);
    this.reference_sequence = referenceFile;
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
    this.nt = 4

    // make sure we have the right intervals
    this.BTIMR = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
  }

  class EvalCalls(@Input vcf: File) extends VariantEval with GATKArgs {
    this.o = swapExt(vcf, ".vcf", ".vcf.eval")
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    this.rodBind :+= RodBind("compOMNI", "VCF", OMNI_VCF)
    this.rodBind :+= RodBind("compHapMap3", "VCF", HM3_VCF)
  }

  class AnnotateCalls(@Input vcf: File, @Output file: File) extends VariantAnnotator with GATKArgs {
    this.o = file
    this.rodBind :+= RodBind("variant", "VCF", vcf)
    this.rodBind :+= RodBind("compOMNI", "VCF", OMNI_VCF)
    this.rodBind :+= RodBind("compHapMap3", "VCF", HM3_VCF)
    this.rodBind :+= RodBind("panel", "VCF", panelVCF)
    this.expression = List("panel.AC", "panel.AN", "panel.AF")
  }

  class MakeTable(@Input vcf: File) extends VariantsToTable with GATKArgs {
    @Output val table = new File(swapExt(vcf, ".vcf", ".vcf.table"))
    this.o = table
    this.rodBind :+= RodBind("variants", "VCF", vcf)
    this.allowMissingData = true
    this.fields = List("CHROM", "POS", "REF", "ALT", "TRANSITION", "HapMap3",
      "OMNI", "AC", "AN", "AF", "panel.AC", "panel.AN", "panel.AF")
  }

  def script = {
    val genotyped = new File(swapExt(bam, ".bam", "_genotyped_at." + panelVCF.getName).getName)
    val panelAnnotated = new File(swapExt(panelVCF, ".vcf", ".annotated.vcf"))
    val annotated = new File(swapExt(genotyped, ".vcf", ".annotated.vcf"))

    add(new GenotypeBAMAtSites(bam, panelVCF, genotyped))
    add(new AnnotateCalls(panelVCF, panelAnnotated))
    add(new AnnotateCalls(genotyped, annotated))
    add(new EvalCalls(annotated))
    add(new MakeTable(annotated))
    add(new MakeTable(panelAnnotated))
  }
}