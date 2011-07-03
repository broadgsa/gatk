package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.RodBind
import org.broadinstitute.sting.queue.extensions.gatk._

class StandardVariantEvaluation extends QScript {
  // todo -- update to released version when things stabilize
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "intervals", doc="intervals to evaluate.  Only supports evaluation on chromosome 20 now, as most evaluation data is there", required=false)
  val TARGET_INTERVAL: String = "20"

  @Argument(shortName = "includeUnion", doc="If provided, we'll create a union of the evaluation data sets for evaluation", required=false)
  val CREATE_UNION: Boolean = false

  @Argument(shortName = "dataDir", doc="Path to the standard evaluation data files", required=false)
  val DATA_DIR = "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/"

  @Argument(shortName = "evalStandard1000GCalls", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val EVAL_STANDARD_1000G_CALLS: Boolean = false

  val COMPS_DIR = DATA_DIR + "/comps/"
  val EVALS_DIR = DATA_DIR + "/evals/"

  @Argument(shortName = "moreSNPsToEval", doc="Path to additional SNP call sets for evaluation", required=false)
  val moreSNPsToEval: List[File] = Nil

  @Argument(shortName = "moreIndelsToEval", doc="Path to additional Indel call sets for evaluation", required=false)
  val moreIndelsToEval: List[File] = Nil

  val VARIANT_TYPES: List[String] = List("indels", "snps")
  val VARIANT_TYPE_VT: Map[String, List[org.broad.tribble.util.variantcontext.VariantContext.Type]] = Map(
    "indels" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.INDEL, org.broad.tribble.util.variantcontext.VariantContext.Type.MIXED, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION),
    "snps" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.SNP, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION)
  )

  val SITES_DIR: String = "sitesFiles"

  // path to b37 DBSNP
  @Argument(shortName = "dbsnp", doc="Path to DBSNP **VCF** for evaluation", required=false)
  val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  //val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf");

  class Comp(val name: String, val evalType: String, val filename: String, val MakeHomVar: Boolean = false) {
    val originalFile = new File(COMPS_DIR + filename)
    val file: File = if ( MakeHomVar ) swapExt(originalFile, ".vcf",".homvar.vcf") else originalFile
    val sitesFile = new File(SITES_DIR + "/" + swapExt(file, ".vcf", ".sites.vcf").getName)
  }

  class Eval(val name: String, val evalType: String, val filename: String, val overrideFile: File = null ) {
    val file: File = if ( overrideFile != null ) overrideFile else new File(EVALS_DIR + "/" + filename)
  }

  var COMPS: List[Comp] = Nil
  def addComp(comp: Comp) { COMPS = comp :: COMPS }

  var EVALS: List[Eval] = Nil
  def addEval(eval: Eval) { EVALS = eval :: EVALS }
  def addEvalFromCMD(file: File, t: String) { addEval(new Eval(file.getName, t, null, file)) }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.intervalsString = List(TARGET_INTERVAL);
    this.reference_sequence = referenceFile;
    this.memoryLimit = 2
  }

  def initializeStandardDataFiles() = {
    //
    // Standard evaluation files for indels
    //
    addComp(new Comp("NA12878.homvar.GATK", "indels", "Indels.NA12878_WGS.filtered_Q50.0_QD5.0_SB-1.0_HR18.vcf", true))
    addComp(new Comp("CG.38samples", "indels", "CG.Indels.leftAligned.b37.vcf"))
    addComp(new Comp("NA12878.homvar.CG", "indels", "NA12878.CG.b37.indels.vcf", true))
    addComp(new Comp("g1k.pilot1.validation", "indels", "pilot1_indel_validation_2009.b37.vcf"))
    addComp(new Comp("NA12878.hand_curated", "indels", "NA12878.validated.curated.polymorphic.indels.vcf"))
    addComp(new Comp("NA12878.Mullikin", "indels", "NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"))


    //
    // INDEL call sets
    //
    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("dindel", "indels", "20110208.chr20.dindel2.EUR.sites.vcf"))
      addEval(new Eval("si", "indels", "20101123.chr20.si.v2.EUR.sites.vcf"))
      addEval(new Eval("gatk", "indels", "EUR.phase1.chr20.broad.filtered.indels.sites.vcf"))
    }

    //
    // Standard evaluation files for SNPs
    //
    addComp(new Comp("NA12878.homvar.GATK", "snps", "NA12878.HiSeq19.cut.vcf", true))
    addComp(new Comp("CG.38samples", "snps", "CG.38samples.b37.vcf"))
    addComp(new Comp("NA12878.homvar.CG", "snps", "NA12878.CG.b37.snps.vcf", true))
    addComp(new Comp("HapMap3.3", "snps", "hapmap3.3.sites_r27_nr.b37_fwd.vcf"))
    addComp(new Comp("OMNI.2.5M", "snps", "omni2.5.1212samples.b37.sites.chr20.monoAreAC0.vcf"))
    addComp(new Comp("g1k.pilot1.validation", "snps", "1000G.snp.validation.b37.vcf"))

    //
    // SNP call sets
    //
    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("1000G.gatk.eurPlus.phase1", "snps", "EUR+.phase1.chr20.broad.recal.vrcut1p0.sites.vcf"))
      addEval(new Eval("1000G.high_specificity.phase1", "snps", "ALL.phase1.chr20.projectConsensus.highSpecificity.snps.genotypes.sites.vcf"))
    }
  }

  def script = {
    val sitesDir = new File(SITES_DIR)
    if ( ! sitesDir.exists ) sitesDir.mkdirs()

    initializeStandardDataFiles();

    // add additional files for evaluation, if necessary
    moreSNPsToEval.foreach(addEvalFromCMD(_, "snps"))
    moreIndelsToEval.foreach(addEvalFromCMD(_, "indels"))

    //
    // create hom-var versions of key files
    //
    for ( comp <- COMPS )
      if ( comp.MakeHomVar )
        add(new SelectHomVars(comp.originalFile, comp.file))

    for ( comp <- COMPS )
        add(new JustSites(comp.file, comp.sitesFile))

    //
    // Loop over evaluation types
    //
    for ( evalType <- VARIANT_TYPES ) {
      var evalsOfType = EVALS.filter(_.evalType == evalType)
      val compsOfType = COMPS.filter(_.evalType == evalType)

      if ( evalsOfType.size > 0 ) {

        // if desired and possible, create a union.X.vcf file
        if ( CREATE_UNION && evalsOfType.size > 1 ) {
          val union: File = new File("union.%s.vcf".format(evalType))
          add(new MyCombine(evalsOfType.map(_.file), union));
          evalsOfType = new Eval("union", evalType, null, union) :: evalsOfType
        }

        // our root VE
        val VE = new MyEval()
        VE.VT = VARIANT_TYPE_VT(evalType)
        VE.o = new File(evalType + ".eval")

        // add evals
        for ( calls <- evalsOfType )
          VE.rodBind :+= RodBind("eval_" + calls.name, "VCF", calls.file)

        // add comps
        //VE.rodBind :+= RodBind("dbsnp", "VCF", MY_DBSNP)
        for ( comp <- compsOfType )
          VE.rodBind :+= RodBind("comp_" + comp.name, "VCF", comp.sitesFile)

        add(VE)
      }
    }
  }

  /**
   * Select homozygous non-reference sites from a single deep data set
   */
  class SelectHomVars(@Input(doc="foo") vcf: File, @Output(doc="foo") out: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("variant", "VCF", vcf)
    this.o = out
    this.select ++= List("\"AC == 2\"")
  }

  /**
   * A simple union
   */
  class MyCombine(@Input(doc="foo") vcfs: List[File], @Output(doc="foo") out: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    for ( vcf <- vcfs )
      this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    this.o = out
  }

  /**
   * A command line (cut) that removes all genotyping information from a file
   */
  class JustSites(@Input(doc="foo") in: File, @Output(doc="foo") out: File) extends CommandLineFunction {
    def commandLine = "cut -f 1-8 %s > %s".format(in, out)
  }

  /**
   * Base class for VariantEval used here
   */
  class MyEval() extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.noST = true
    this.evalModule :+= "ValidationReport"
  }
}

