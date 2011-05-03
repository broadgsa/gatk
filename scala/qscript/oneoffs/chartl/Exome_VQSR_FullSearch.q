import collection.mutable.HashMap
import java.io.{PrintWriter, PrintStream}
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, ScatterFunction, GatherFunction}
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.interval.IntervalSetRule
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._

class Exome_VQSR_FullSearch extends QScript {
  qScript =>

  // COMMAND LINE ARGUMENTS

  // VARIABLES USED
  val SCRIPT_BASE_NAME = "Exome_VQSR_Search"
  val UG_CALL_THRESH = 10.0;
  val VQSR_CALL_THRESH = List(10.0,20.0,30.0,40.0,50.0)
  val VQSR_ANNOTATIONS_TO_USE = List(List("QD","SB"),List("QD","SB","HRun"),List("QD","SB","HaplotypeScore"),List("QD","SB","HRun","HaplotypeScore"))
  val HM3_SITES = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf")
  val OMNI_CHIP = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf")
  val AXIOM_CHIP = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/affymetrix_axiom/Affymetrix_Axiom_DB_2010_v4_b37.noOmni.noHM3.vcf")
  val DBSNP_129 = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf")
  val SENSITIVITY = (new Range(0,19,1)).map(u => 90.0 + 0.5*u).toList
  val RECALIBRATE_TOGETHER = List(true,false)
  val VQSR_HAPMAP_PRIOR = "15.0"
  val VQSR_OMNI_PRIOR = "12.0"

  var VQSR_RODBINDS : HashMap[String,List[RodBind]] = new HashMap[String,List[RodBind]]
  val VQSR_TAG_FT = "known=false,training=true,truth=%s,prior=%s"
  val VQSR_DBSNP_TAG = "known=true,training=false,truth=false,prior=0.1"


  for ( tf <- List( (true,false),(false,true),(true,true)) ) {
    var mrb: List[RodBind] = Nil
    val ext = (if ( tf._1 ) "HT" else "HF") + (if ( tf._2 ) "OT" else "OF")
    val hmSt = if ( tf._1 ) "true" else "false"
    val omSt = if ( tf._2 ) "true" else "false"
    mrb :+= RodBind("dbsnp","VCF",DBSNP_129,VQSR_DBSNP_TAG)
    mrb :+= RodBind("HapMap3","VCF",HM3_SITES,VQSR_TAG_FT.format(hmSt,VQSR_HAPMAP_PRIOR))
    mrb :+= RodBind("Omni","VCF",OMNI_CHIP,VQSR_TAG_FT.format(omSt,VQSR_OMNI_PRIOR))
    VQSR_RODBINDS += new Tuple2(ext,mrb)
  }

  val BAM_FILES : List[File] = asScalaIterator((new XReadLines(new File("/humgen/gsa-hphome1/chartl/projects/oneoffs/VQSR_Exome/resources/broad.bam.list")))).map(u => new File(u)).toList
  val REF : File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
  val INTS : File = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
  val EXPAND_INTS = 40
  val HAND_FILTERED : File = new File("/humgen/1kg/exomes/results/broad.wex.96samples/v1/1KGBroadWEx.variants.vcf")
  val GATK_JAR : File = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")

  def script() {
    trait CommandLineGATKArgs extends CommandLineGATK {
      this.intervals :+= INTS
      this.jarFile = GATK_JAR
      this.reference_sequence = REF
      this.memoryLimit = Some(4)
    }

    val ei : ExpandIntervals = new ExpandIntervals(INTS,1,EXPAND_INTS, new File("Resources", SCRIPT_BASE_NAME + ".flanks.interval_list"), REF, "INTERVALS", "INTERVALS")
    ei.jobOutputFile = new File(".queue/logs/Overall/ExpandIntervals.out")

    if (EXPAND_INTS > 0) {
      //add(ei)
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (EXPAND_INTS > 0) {
        this.intervals :+= ei.outList.getAbsoluteFile
      }
    }

    val callSNPsAndIndels = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    callSNPsAndIndels.analysisName = SCRIPT_BASE_NAME+"_calls"
    callSNPsAndIndels.jobOutputFile = new File(".queue/logs/SNPCalling/UnifiedGenotyper.snps.out")
    callSNPsAndIndels.memoryLimit = Some(6)
    callSNPsAndIndels.downsample_to_coverage = Some(600)
    callSNPsAndIndels.input_file = BAM_FILES
    callSNPsAndIndels.rodBind :+= RodBind("dbsnp", "vcf", DBSNP_129)
    callSNPsAndIndels.out = new File(SCRIPT_BASE_NAME+".rawCalls.vcf")
    callSNPsAndIndels.stand_call_conf = Some(UG_CALL_THRESH)

    callSNPsAndIndels.scatterCount = 50
    callSNPsAndIndels.setupScatterFunction = {
      case scatter: ScatterFunction =>
        scatter.commandDirectory = new File("SnpCalls/ScatterGather")
        scatter.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Scatter.out")
    }
    callSNPsAndIndels.setupCloneFunction = {
      case (clone: CloneFunction, index: Int) =>
        clone.commandDirectory = new File("SnpCalls/ScatterGather/Scatter_%s".format(index))
        clone.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Scatter_%s.out".format(index))
    }
    callSNPsAndIndels.setupGatherFunction = {
      case (gather: GatherFunction, source: ArgumentSource) =>
        gather.commandDirectory = new File("SnpCalls/ScatterGather/Gather_%s".format(source.field.getName))
        gather.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Gather_%s.out".format(source.field.getName))
    }

    //add(callSNPsAndIndels)

    class ExtractSNPs(in: File, out: File) extends InProcessFunction {

      @Input(doc="foo")
      var inputVCF : File = in
      @Output(doc="foo")
      var outputVCF : File = out

      def canPrint(line: String) : Boolean = {
        //System.out.println(line)
        if ( line.startsWith("#") ) { return true }
        val spline = line.split("\t",6)
        //System.out.println(spline.apply(3)+"  "+spline.apply(4))
        if ( spline.apply(3).size > 1 || spline.apply(4).size > 1 ) { return false }
        true
      }

      def run() {
        val outWriter = new PrintWriter(new PrintStream(outputVCF))
        asScalaIterator(new XReadLines(inputVCF)).foreach(u => {
          if ( canPrint(u) ) {
            outWriter.println(u)
          }
        })

        outWriter.close
      }
    }

    val extractSNPs : ExtractSNPs = new ExtractSNPs(callSNPsAndIndels.out,new File(SCRIPT_BASE_NAME+".snpCalls.vcf"))
    //add(extractSNPs)

    def getPath(annoList: List[String], jointRecal: Boolean) : File = {
      new File("VQSR/%s/%s".format( annoList.reduceLeft( _ + "." + _ ) , if(jointRecal) "together" else "separate")  )
    }

    var filterMap : HashMap[Double,File] = new HashMap[Double,File]

    for ( thresh <- VQSR_CALL_THRESH ) {
      var filterQual = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
      filterQual.rodBind :+= new RodBind("variant","VCF",extractSNPs.outputVCF.getAbsoluteFile)
      filterQual.filterExpression :+= "\'QUAL < %.1f\'".format(thresh)
      filterQual.filterName :+= "LowQual"
      filterQual.out = new File(SCRIPT_BASE_NAME+".filterQual%.1f.vcf".format(thresh))
      add(filterQual)
      filterMap += new Tuple2(thresh,filterQual.out.getAbsoluteFile)
    }

    for ( annotations <- VQSR_ANNOTATIONS_TO_USE ) {
      for ( recalTogether <- RECALIBRATE_TOGETHER ) {
        val directory = getPath(annotations,recalTogether)
        for ( call_thresh <- VQSR_CALL_THRESH ) {
          for ( vqsr_rb <- VQSR_RODBINDS.iterator ) {
            trait VQSR_Args extends VariantRecalibrator {
              this.allPoly = true
              this.analysisName = "VQSR_%s_%s_%.1f".format( annotations.reduceLeft( _ + "." + _), if ( recalTogether ) "true" else "false", call_thresh)
              this.commandDirectory = directory
              this.use_annotation ++= annotations
              this.tranche ++= SENSITIVITY.map(u => "%.1f".format(u))
              this.rodBind :+= RodBind("inputData","VCF",filterMap.get(call_thresh).get)
              this.rodBind ++= vqsr_rb._2
              this.memoryLimit = Some(8)
            }
            val nameFormat = SCRIPT_BASE_NAME+".%1f.%s.".format(call_thresh,vqsr_rb._1)+"%s."
            if ( recalTogether ) {
              var vqsr = new VariantRecalibrator with VQSR_Args with ExpandedIntervals with CommandLineGATKArgs
              vqsr.tranchesFile = new File(nameFormat.format("both")+"tranche")
              vqsr.recalFile = new File(nameFormat.format("both")+"recal")
              add(vqsr)
              addAll(eval(vqsr, ei.outList, "flanks"))
              addAll(eval(vqsr, INTS, "exons"))
            } else {
              var exons = new VariantRecalibrator with VQSR_Args with CommandLineGATKArgs
              exons.tranchesFile = new File(nameFormat.format("exons")+"tranche")
              exons.recalFile = new File(nameFormat.format("exons")+"recal")
              var flanks = new VariantRecalibrator with VQSR_Args
              flanks.intervals :+= ei.outList.getAbsoluteFile
              flanks.jarFile = GATK_JAR
              flanks.memoryLimit = Some(8)
              flanks.reference_sequence = REF
              flanks.tranchesFile = new File(nameFormat.format("flanks")+"tranche")
              flanks.recalFile = new File(nameFormat.format("flanks")+"recal")
              add(exons,flanks)
              addAll(eval(exons))
              addAll(eval(flanks))
            }
          }
        }
      }
    }
  }

  // want to apply and eval
  def eval(recal: VariantRecalibrator) : List[QFunction] = { eval(recal,null,"") }
  def eval(recal: VariantRecalibrator, list: File, ext: String) : List[QFunction] = {
    var functions : List[QFunction] = Nil
    trait ImplicitArgs extends CommandLineGATK {
      this.jarFile = recal.jarFile
      this.reference_sequence = recal.reference_sequence
      this.commandDirectory = recal.commandDirectory
      if ( list == null ) {
        this.intervals ++= recal.intervals
      } else {
        this.intervals :+= list
      }
    }

    trait ApplyArgs extends ApplyRecalibration with ImplicitArgs {
      this.tranchesFile = recal.tranchesFile
      this.recalFile = recal.recalFile
      for ( r <- recal.rodBind ) {
        if ( r.trackName.startsWith("input") ) {
          this.rodBind :+= r
        }
      }
      this.memoryLimit = Some(4)
    }

    trait EvalArgs extends VariantEval with ImplicitArgs {
      this.stratificationModule = List("Novelty")
      this.evalModule = List("TiTvVariantEvaluator","CountVariants","GenotypeConcordance")
      this.rodBind :+= RodBind("dbsnp","VCF",DBSNP_129)
      this.rodBind :+= RodBind("compAxiom","VCF",AXIOM_CHIP)
      this.memoryLimit = Some(4)
    }

    val extender = if ( ext != null ) ".cut%.1f."+ext else ".cut%.1f"
    for ( sens <- SENSITIVITY ) {
      var cut = new ApplyRecalibration with ApplyArgs
      cut.analysisName = recal.analysisName+extender.format(sens)
      val vcfExt = extender.format(sens)+".vcf"
      cut.out = swapExt(cut.recalFile,".recal",vcfExt)
      cut.ts_filter_level = Some(sens)
      functions :+= cut

      var eval = new VariantEval with EvalArgs
      eval.analysisName = cut.analysisName+".eval"
      eval.out = swapExt(cut.out,".vcf",".eval")
      eval.rodBind :+= RodBind("evalContrastive","VCF",cut.out)
      functions :+= eval
    }

    functions
  }
}
