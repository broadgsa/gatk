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

  @Argument(shortName = "truth", doc="VQSR truth file", required=false)
  var truthFile: File = new File("/humgen/gsa-scr1/delangel/devine_data/indel_hg19_051711_leftAligned_75percent_chr20.vcf"  )
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
    this.jobQueue = "week"
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
    addComp(new Comp("Phase1Validation", "indels", "/humgen/gsa-scr1/delangel/VQSRIndels/1KG_Validation_Phase1_SNPs_05032011.HG19.finalized.vcf"))


    //
    // INDEL call sets
    //

    if ( EVAL_STANDARD_1000G_CALLS ) {
//      addEval(new Eval("dindel", "indels", "/humgen/gsa-scr1/delangel/officialCalls/20110201_chr20_phase1_indels/dindel/20110208.chr20.dindel2.ALL.sites.vcf"))
      addEval(new Eval("si", "indels", "/humgen/gsa-scr1/delangel/officialCalls/20101123.chr20.si.v2.combined.sites.leftAligned.vcf"))
      addEval(new Eval("bi", "indels", "/humgen/1kg/processing/official_release/phase1/ALL.wgs.broad.20101123.indels.sites.vcf"))
      addEval(new Eval("bc", "indels", "/humgen/gsa-scr1/delangel/officialCalls/20110201_chr20_phase1_indels/ALL.chr20.bc.20101123.indels.sites.leftAligned.vcf"))
      addEval(new Eval("ox", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.chr20.Oxford.20110407.indels.genotypes.sites.vcf"))
      addEval(new Eval("2of5", "indels", "/humgen/gsa-scr1/delangel/otherIndelCallerAnalysis/ALL.indels.2of5.chr20.vcf"))
      }

    //
    // Standard evaluation files for SNPs
    //
    /*   addComp(new Comp("NA12878.homvar.GATK", "snps", "NA12878.HiSeq19.cut.vcf", true))
    addComp(new Comp("CG.38samples", "snps", "CG.38samples.b37.vcf"))
    addComp(new Comp("NA12878.homvar.CG", "snps", "NA12878.CG.b37.snps.vcf", true))
    addComp(new Comp("HapMap3.3", "snps", "hapmap3.3.sites_r27_nr.b37_fwd.vcf"))
    addComp(new Comp("OMNI.2.5M", "snps", "omni2.5.1212samples.b37.sites.chr20.monoAreAC0.vcf"))
    addComp(new Comp("g1k.pilot1.validation", "snps", "1000G.snp.validation.b37.vcf"))
    */
    //
    // SNP call sets
    //
  }

  def script = {

    initializeStandardDataFiles();

    // add additional files for evaluation, if necessary
    //moreSNPsToEval.foreach(addEvalFromCMD(_, "snps"))
    //moreIndelsToEval.foreach(addEvalFromCMD(_, "indels"))


    var ts:Double = 0.0
    var tranches =  List("99.9","99.0","98.0","97.0","95.0","92.0","90.0")

    var numG:Int = qscript.numG
    var pctBad:Double = qscript.pctBad
    val runName:String = qscript.runName +  "_mG%d_pb%1.2f_QD_FS_HS_RP_IC".format(numG,pctBad)

    for( pop <- qscript.populations ) {
      val rawCalls = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110516Dev/calls/chr20/%s/%s.phase1.chr20.raw.indels.vcf".format(pop,pop))
      //val filteredCalls = new File(baseName + ".v5."+runStr+".filtered.vcf")
      // val runname =
      //val clusterFile = new File(baseName + ".omni.clusters")
      //val recalibratedCalls = new File(baseName + ".recal.vcf")
      var tranchesFile = new File(qscript.baseDir +"%s_%s.tranches".format(pop,runName))
      var recalFile = new File(qscript.baseDir +"%s_%s.recal".format(pop,runName))
      var rscriptFile = new File(qscript.baseDir +"%s_%s.plots.R".format(pop,runName))



      var vr = new VariantRecalibrator with CommandLineGATKArgs
      vr.rodBind :+= RodBind("input", "VCF",rawCalls )
      vr.rodBind :+= RodBind("truth", "VCF",qscript.truthFile,"known=true,training=true,truth=true,prior=20.0" )
      vr.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
      vr.tranchesFile = tranchesFile
      vr.recalFile = recalFile
      vr.rscriptFile = rscriptFile
      vr.an = List("QD","FS","HaplotypeScore","ReadPosRankSum","InbreedingCoeff")
      vr.maxGaussians = Some(numG)
      vr.tranche = tranches
      vr.nt = Some(8)
      vr.percentBad = Some(pctBad)
      add(vr)


      for (tas: String <- tranches) {
        ts = tas.toDouble
        val outFile = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110516Dev/calls/chr20/%s/%s.phase1.chr20.recal_%s_ts_%4.1f.indels.sites.vcf".format(pop,pop,runName,ts))

        var ar = new ApplyRecalibration with CommandLineGATKArgs
        ar.rodBind :+= RodBind("input", "VCF",rawCalls )
        ar.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
        ar.tranchesFile = tranchesFile
        ar.recalFile = recalFile
        ar.ts_filter_level = Some(ts)
        ar.sites_only = true
        ar.o = outFile
        add(ar)
      }
    }


    val VE = new MyEval()
    VE.VT = VARIANT_TYPE_VT("indels")
    VE.o = new File(OUT_DIR+"/"+ runName + ".eval")
    //VE.nt = Some(8)

    for (tas: String <- tranches) {
      ts = tas.toDouble
      var cm = new CombineVariants with CommandLineGATKArgs

      cm.o = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110516Dev/calls/chr20/ALL.phase1.chr20.recal_%s_ts_%4.1f.indels.sites.vcf".format(runName,ts))
      for( pop <- qscript.populations ) {
        val outFile = new File("/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110516Dev/calls/chr20/%s/%s.phase1.chr20.recal_%s_ts_%4.1f.indels.sites.vcf".format(pop,pop,runName,ts))
        cm.rodBind :+= RodBind(pop, "VCF", outFile)
      }
      add(cm)

      VE.rodBind :+= RodBind("eval_ts%4.1f".format(ts), "VCF", cm.o)

    }
    // add evals
    for ( calls <- EVALS )
      VE.rodBind :+= RodBind("eval_" + calls.name, "VCF", calls.file)

    // add comps
//    VE.rodBind :+= RodBind("dbsnp", "VCF", MY_DBSNP)
    for ( comp <- COMPS )
      VE.rodBind :+= RodBind("comp_" + comp.name, "VCF", comp.file)

    add(VE)

   }


  /**
   * Base class for VariantEval used here
   */
  class MyEval() extends VariantEval with CommandLineGATKArgs {
    this.noST = true
    this.evalModule :+= "ValidationReport"
    //this.evalModule :+= "IndelMetricsByAC"
    this.evalModule :+= "IndelStatistics"
    this.evalModule :+= "CountVariants"
    //this.evalModule :+= "IndelClasses"
  }



}
