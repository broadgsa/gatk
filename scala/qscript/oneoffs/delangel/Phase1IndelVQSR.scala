import collection.SeqLike._
import management.CompilationMXBean
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.IndelStatistics
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.RodBind._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils
import scala.Some
;

class Phase1IndelVQSR extends QScript {
  qscript =>
  // todo -- update to released version when things stabilize
  @Argument(shortName = "gatk",doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/humgen/gsa-scr1/delangel/Sting_dev/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "intervals", doc="intervals to evaluate.  Only supports evaluation on chromosome 20 now, as most evaluation data is there", required=false)
  val TARGET_INTERVAL: String = "20"

  @Argument(shortName = "dataDir", doc="Path to the standard evaluation data files", required=false)
  val DATA_DIR = "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/"
  @Argument(shortName = "baseDir", doc="Path to the standard evaluation data files", required=false)
  val baseDir = "/humgen/gsa-scr1/delangel/VQSRIndels/data/"

  @Argument(shortName = "outDir", doc="Path to the output files", required=false)
  val OUT_DIR = "/humgen/gsa-scr1/delangel/VQSRIndels"

  @Argument(shortName = "rawCalls", doc="VQSR raw input file", required=false)
  var rawCalls: File = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110525VQSRConsensus/calls/combined.phase1.chr20.raw.indels.vcf")

  @Argument(shortName = "truth", doc="VQSR truth file", required=false)
  var truthFile: File = new File("/humgen/gsa-scr1/delangel/devine_data/indel_hg19_051711_leftAligned_75percent_chr20.vcf"  )

  @Argument(shortName = "training", doc="VQSR training file", required=false)
  var trainingFile: File = new File("/humgen/gsa-scr1/delangel/devine_data/indel_hg19_051711_leftAligned_75percent_chr20.vcf"  )

  var noMultiallelicSites: Boolean = false;

  val populations = List("EUR","AMR","ASN","AFR")

  @Argument(shortName = "evalStandard1000GCalls", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val EVAL_STANDARD_1000G_CALLS: Boolean = true

  @Argument(shortName = "numG", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val numG: Int = 4

  @Argument(shortName = "pctBad", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val pctBad: Double = 0.05

  @Argument(shortName = "runName", doc="Run Name", required=false)
  val runName:String = "mills100"
  val COMPS_DIR = DATA_DIR + "comps/"
  val EVALS_DIR = DATA_DIR + "evals/"

  @Argument(shortName = "createAllPos", doc="If provided, create all POPS file", required=false)
  val CREATE_ALL_POPS_FILE: Boolean = false

  @Argument(shortName = "pops", doc="Populations to do", required=false)
  val moreIndelsToEval: List[String] = List("EUR","ASN","AFR","AMR")


  val VARIANT_TYPES: List[String] = List("indels", "snps")

  val VARIANT_TYPE_VT: Map[String, List[org.broad.tribble.util.variantcontext.VariantContext.Type]] = Map(
    "indels" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.INDEL, org.broad.tribble.util.variantcontext.VariantContext.Type.MIXED, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION),
    "snps" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.SNP, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION)
  )

  val SITES_DIR: String = "sitesFiles"

  // path to b37 DBSNP
  val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf"


  val DINDEL: String =  "/humgen/gsa-scr1/delangel/officialCalls/20110201_chr20_phase1_indels/dindel/20110208.chr20.dindel2.ALL.sites.fixed.vcf"
  val SI: String =  "/humgen/gsa-scr1/delangel/officialCalls/20101123.chr20.si.v2.combined.sites.leftAligned.vcf"
  val BI: String =   "/humgen/1kg/processing/official_release/phase1/ALL.wgs.broad.20101123.indels.sites.vcf"
  val BC: String =   "/humgen/gsa-scr1/delangel/officialCalls/20110201_chr20_phase1_indels/ALL.chr20.bc.20101123.indels.sites.leftAligned.vcf"
  val OX: String =  "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.chr20.Oxford.20110407.indels.genotypes.sites.vcf"

  var COMPS: List[Comp] = Nil
  def addComp(comp: Comp) { COMPS = comp :: COMPS }

  var EVALS: List[Eval] = Nil
  def addEval(eval: Eval) { EVALS = eval :: EVALS }
  def addEvalFromCMD(file: File, t: String) { addEval(new Eval(file.getName, t, file.getName)) }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = Some(2)
    // this.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )
    this.jobQueue = "hour"
    this.intervalsString = List(TARGET_INTERVAL);

  }
  class Comp(val name: String, val evalType: String, val filename: String) {
    val file: File = new File(filename)
  }

  class Eval(val name: String, val evalType: String, val filename: String ) {
    val file: File = new File(filename)
  }

  def initializeStandardDataFiles() = {
    //
    // Standard evaluation files for indels
    //
    //addComp(new Comp("CG.38samples", "indels", COMPS_DIR+"CG.Indels.leftAligned.b37.vcf"))
    addComp(new Comp("g1k.pilot1.validation", "indels", COMPS_DIR+"pilot1_indel_validation_2009.b37.vcf"))
    //addComp(new Comp("NA12878.hand_curated", "indels", "NA12878.validated.curated.polymorphic.indels.vcf"))
    addComp(new Comp("NA12878.Mullikin", "indels", COMPS_DIR+"NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"))
    addComp(new Comp("Mills.25pct", "indels", "/humgen/gsa-scr1/delangel/devine_data/indel_hg19_051711_leftAligned_25percent_chr20.vcf"))
    //addComp(new Comp("Phase1Validation", "indels", "/humgen/gsa-scr1/delangel/VQSRIndels/1KG_Validation_Phase1_SNPs_05032011.HG19.finalized.vcf"))
   addComp(new Comp("Phase1Validation", "indels", "/humgen/gsa-scr1/delangel/VQSRIndels/1000G.20101123.validation_set_v1.QCed.indels.vcf"))


    //
    // INDEL call sets
    //

    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("dindel", "indels",qscript.DINDEL))
      addEval(new Eval("si", "indels",qscript.SI))
      addEval(new Eval("bi", "indels", qscript.BI))
      addEval(new Eval("bc", "indels", qscript.BC))
      addEval(new Eval("ox", "indels", qscript.OX))
      addEval(new Eval("2of5", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.2of5.chr20.vcf"))
      addEval(new Eval("2of5noMulti", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.2of5.chr20.noMultiAllelic.vcf"))
      addEval(new Eval("union", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.combined.chr20.vcf"))
//       addEval(new Eval("unionNoMulti", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.combined.chr20.noMultiAllelic.vcf"))
    }

  }

  def script = {

    initializeStandardDataFiles();

    var ts:Double = 0.0
    var tranches =  List("99.9","99.0","98.0","97.0","96.0","95.0","92.0","90.0","85.0","80.0","70.0")

    var numG:Int = qscript.numG
    var pctBad:Double = qscript.pctBad
    val runName:String = qscript.runName +  "_mG%d_pb%1.2f_QD_FS_HS_RP_IC".format(numG,pctBad)

    val rawCalls = qscript.rawCalls
    var tranchesFile = new File(qscript.baseDir +"%s.tranches".format(runName))
    var recalFile = new File(qscript.baseDir +"%s.recal".format(runName))
    var rscriptFile = new File(qscript.baseDir +"%s.plots.R".format(runName))



    var vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.rodBind :+= RodBind("input", "VCF",rawCalls )
    vr.rodBind :+= RodBind("truth", "VCF",qscript.truthFile,"known=true,training=false,truth=true,prior=15.0" )
    vr.rodBind :+= RodBind("training", "VCF",qscript.trainingFile,"known=true,training=true,truth=false,prior=12.0" )
    //vr.rodBind :+= RodBind("training2", "VCF",qscript.dindelCalls,"known=true,training=true,truth=false,prior=12.0" )
    vr.rodBind :+= RodBind("dbsnp", "VCF",qscript.MY_DBSNP,"known=true,training=false,truth=false,prior=8.0" )

    vr.rodBind :+= RodBind("BC", "VCF",qscript.BC,"consensus=true" )
    vr.rodBind :+= RodBind("BI", "VCF",qscript.BI,"consensus=true" )
    vr.rodBind :+= RodBind("SI", "VCF",qscript.SI,"consensus=true" )
    vr.rodBind :+= RodBind("DINDEL", "VCF",qscript.DINDEL,"consensus=true" )
    vr.rodBind :+= RodBind("OXFORD", "VCF",qscript.OX,"consensus=true" )

    vr.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
    vr.tranchesFile = tranchesFile
    vr.recalFile = recalFile
    vr.rscriptFile = rscriptFile
//    vr.an = List("QD","FS","HaplotypeScore","ReadPosRankSum","InbreedingCoeff","SB","")
    vr.an = List("QD","FS","HaplotypeScore","ReadPosRankSum","InbreedingCoeff")
    vr.maxGaussians = Some(numG)
    vr.tranche = tranches
    vr.nt = Some(8)
    vr.percentBad = Some(pctBad)
    vr.std = Some(12.0)
    //vr.ignore_filter = List("LowQual")
    add(vr)

    val VE = new MyEval()
    VE.VT = VARIANT_TYPE_VT("indels")
    VE.o = new File(OUT_DIR+"/"+ runName + ".eval")

    for (tas: String <- tranches) {
      ts = tas.toDouble
      val outFile = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110603VQSRConsensus/calls/phase1.chr20.recal_%s_ts_%4.1f.indels.sites.vcf".format(runName,ts))

      var ar = new ApplyRecalibration with CommandLineGATKArgs
      ar.rodBind :+= RodBind("input", "VCF",rawCalls )
      ar.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
      ar.tranchesFile = tranchesFile
      ar.recalFile = recalFile
      ar.ts_filter_level = Some(ts)
      ar.sites_only = true
      ar.o = outFile
      add(ar)

      VE.rodBind :+= RodBind("eval_ts%4.1f".format(ts), "VCF", ar.o)

    }


    //VE.nt = Some(8)

    // add evals
    for ( calls <- EVALS )
      VE.rodBind :+= RodBind("eval_" + calls.name, "VCF", calls.file)

    // add comps
    //    VE.rodBind :+= RodBind("dbsnp", "VCF", MY_DBSNP)
    for ( comp <- COMPS )
      VE.rodBind :+= RodBind("comp_" + comp.name, "VCF", comp.file)

    add(VE)

    var  ve2 = new MyEval
    for (tas: String <- tranches) {
      ts = tas.toDouble
      val outFile = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110603VQSRConsensus/calls/phase1.chr20.recal_%s_ts_%4.1f.indels.sites.vcf".format(runName,ts))
      ve2.rodBind :+= RodBind("eval_ts%4.1f".format(ts), "VCF", outFile)
    }

    // comps are now other callsets to measure overlap
    ve2.rodBind :+= RodBind("comp_dindel", "VCF",qscript.DINDEL)
    ve2.rodBind :+= RodBind("comp_bc", "VCF", qscript.BC)
    ve2.rodBind :+= RodBind("comp_bi", "VCF", qscript.BI)
    ve2.rodBind :+= RodBind("comp_ox", "VCF", qscript.OX)
    ve2.rodBind :+= RodBind("comp_2of5", "VCF", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.2of5.chr20.vcf")
    ve2.VT = VARIANT_TYPE_VT("indels")
    ve2.o = new File(OUT_DIR+"/"+ runName + ".comps.eval")
    add(ve2)
  }


  /**
   * Base class for VariantEval used here
   */
  class MyEval() extends VariantEval with CommandLineGATKArgs {
    this.noST = true
    this.nt = Some(8)
    this.evalModule :+= "ValidationReport"
    //this.evalModule :+= "IndelMetricsByAC"
    this.evalModule :+= "IndelStatistics"
    this.evalModule :+= "CountVariants"
    this.evalModule :+= "CompOverlap"
    //this.evalModule :+= "IndelClasses"
  }



}
