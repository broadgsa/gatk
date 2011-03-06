package oneoffs.depristo

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.extensions.gatk.RodBind
import org.broadinstitute.sting.queue.extensions.gatk._

class StandardVariantEvalation extends QScript {
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val TARGET_INTERVAL: String = "20";

  val DATA_DIR = "data/"     // todo -- make into an argument

  val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf");

  val JOB_QUEUE = "hour"

  // todo -- add arguments to include other files for SNPs and indels

  val VARIANT_TYPES: List[String] = List("indels", "snps")
  val VARIANT_TYPE_VT: Map[String, List[org.broad.tribble.util.variantcontext.VariantContext.Type]] = Map(
    "indels" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.INDEL, org.broad.tribble.util.variantcontext.VariantContext.Type.MIXED, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION),
    "snps" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.SNP, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION)
  )

  class Comp(val name: String, val evalType: String, val filename: String, val MakeHomVar: Boolean = false) {
    val originalFile = new File(DATA_DIR + filename)
    val file: File = if ( MakeHomVar ) swapExt(originalFile, ".vcf",".homvar.vcf") else originalFile
    val sitesFile = swapExt(file, ".vcf", ".sites.vcf")
  }

  class Eval(val name: String, val evalType: String, val file: File) {}

  var COMPS: List[Comp] = Nil
  def addComp(comp: Comp) { COMPS = comp :: COMPS }

  var EVALS: List[Eval] = Nil
  def addEval(eval: Eval) { EVALS = eval :: EVALS }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.intervalsString = List(TARGET_INTERVAL);
    this.reference_sequence = referenceFile;
    this.jobQueue = JOB_QUEUE;
    this.memoryLimit = Some(2)
  }

  def script = {
    //
    // Standard evaluation files for indels
    //
    addComp(new Comp("homVarNA12878GATK", "indels", "Indels.NA12878_WGS.filtered_Q50.0_QD5.0_SB-1.0_HR18.vcf", true))
    addComp(new Comp("CG38", "indels", "CG.Indels.leftAligned.b37.vcf"))
    addComp(new Comp("homVarNA12878CG", "indels", "NA12878.CG.b37.indels.vcf", true))
    addComp(new Comp("Pilot1Validation", "indels", "pilot1_indel_validation_2009.b37.vcf"))
    addComp(new Comp("Pilot1Validation", "indels", "NA12878.validated.curated.polymorphic.indels.vcf"))

    //
    // INDEL call sets
    //
    addEval(new Eval("dindel", "indels", new File("20110208.chr20.dindel2.EUR.sites.vcf")))
    addEval(new Eval("si", "indels", new File("20101123.chr20.si.v2.EUR.sites.vcf")))
    addEval(new Eval("gatk", "indels", new File("EUR.phase1.chr20.broad.filtered.indels.sites.vcf")))

    //
    // Standard evaluation files for SNPs
    //
    addComp(new Comp("homVarNA12878GATK", "snps", "NA12878.HiSeq19.cut.vcf", true))
    addComp(new Comp("CG38", "snps", "CG.38samples.b37.vcf"))
    addComp(new Comp("homVarNA12878CG", "snps", "NA12878.CG.b37.snps.vcf", true))

    //
    // SNP call sets
    //
    addEval(new Eval("gatk", "snps", new File("EUR+.phase1.chr20.broad.recal.vrcut1p0.sites.vcf")))

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

      if ( evalsOfType.size > 1 ) {
        val union: File = new File("union.%s.vcf".format(evalType))
        add(new MyCombine(evalsOfType.map(_.file), union));
        evalsOfType = new Eval("union", evalType, union) :: evalsOfType
      }

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


  class SelectHomVars(@Input(doc="foo") vcf: File, @Output(doc="foo") out: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("variant", "VCF", vcf)
    this.o = out
    this.select ++= List("\"AC == 2\"")
  }

  class MyCombine(@Input(doc="foo") vcfs: List[File], @Output(doc="foo") out: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    for ( vcf <- vcfs )
      this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    this.o = out
  }

  class JustSites(@Input(doc="foo") in: File, @Output(doc="foo") out: File) extends CommandLineFunction {
    def commandLine = "cut -f 1-8 %s > %s".format(in, out)
    this.jobQueue = JOB_QUEUE;
  }

  class MyEval() extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.noST = true
    this.evalModule :+= "ValidationReport"
  }
}

