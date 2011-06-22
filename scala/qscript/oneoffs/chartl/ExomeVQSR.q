import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines

class FullCallingPipeline extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="G")
  var gatkJar: File = _

  @Input(doc="level of parallelism for UnifiedGenotyper (both for SNPs and indels).  By default is set to 20.", shortName="varScatter", required=false)
  var num_var_scatter_jobs = 20

  @Argument(doc="expand each target in input intervals by the specified number of bases (50 bases by default)", shortName="expand", required=false)
  var expandIntervals = 50

  private var pipeline: Pipeline = _

  private final val picardFixMatesClass = "net.sf.picard.sam.FixMateInformation"


  val BAM_FILES : List[File] = (new XReadLines(new File("/humgen/gsa-hphome1/chartl/projects/oneoffs/VQSR_Exome/resources/broad.bam.list"))).readLines.map(u => new File(u)).toList
  val DBSNP : File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf")
  val REF : File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
  val INTS : File = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
  val BASE : String = "exon_vqsr"

  val handFiltered : File = new File("/humgen/1kg/exomes/results/broad.wex.96samples/v1/1KGBroadWEx.variants.vcf")


  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals :+= INTS
    this.jarFile = qscript.gatkJar
    this.reference_sequence = REF
    this.memoryLimit = Some(4)
  }
  // ------------ SETUP THE PIPELINE ----------- //


  def script = {

    endToEnd(BASE,"cleaned")
  }

  def endToEnd(base: String, bamType: String) = {
    val bamFiles = BAM_FILES

    val ei : ExpandIntervals = new ExpandIntervals(INTS,1,qscript.expandIntervals, new File("Resources", base + ".flanks.interval_list"), REF, "INTERVALS", "INTERVALS")
    ei.jobOutputFile = new File(".queue/logs/Overall/ExpandIntervals.out")

    if (qscript.expandIntervals > 0) {
      //add(ei)
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0) {
        this.intervals :+= ei.outList
      }
    }

    // Call indels
    val indels = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    indels.analysisName = base + "_indels"
    indels.jobOutputFile = new File(".queue/logs/IndelCalling/UnifiedGenotyper.indels.out")
    indels.memoryLimit = Some(6)
    indels.downsample_to_coverage = Some(600)
    indels.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.DINDEL
    indels.input_file = bamFiles
    indels.rodBind :+= RodBind("dbsnp", "vcf", DBSNP)
    indels.out = new File("IndelCalls", base+".indels.vcf")

    indels.scatterCount = qscript.num_var_scatter_jobs
    indels.setupScatterFunction = {
      case scatter: ScatterFunction =>
        scatter.commandDirectory = new File("IndelCalls/ScatterGather")
        scatter.jobOutputFile = new File(".queue/logs/IndelCalling/ScatterGather/Scatter.out")
    }
    indels.setupCloneFunction = {
      case (clone: CloneFunction, index: Int) =>
        clone.commandDirectory = new File("IndelCalls/ScatterGather/Scatter_%s".format(index))
        clone.jobOutputFile = new File(".queue/logs/IndelCalling/ScatterGather/Scatter_%s.out".format(index))
    }
    indels.setupGatherFunction = {
      case (gather: GatherFunction, source: ArgumentSource) =>
        gather.commandDirectory = new File("IndelCalls/ScatterGather/Gather_%s".format(source.field.getName))
        gather.jobOutputFile = new File(".queue/logs/IndelCalling/ScatterGather/Gather_%s.out".format(source.field.getName))
    }

    // Filter indels
    val filteredIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filteredIndels.analysisName = base + "_filteredIndels"
    filteredIndels.jobOutputFile = new File(".queue/logs/IndelCalling/VariantFiltration.indels.out")
    filteredIndels.filterName ++= List("IndelQUALFilter","IndelSBFilter","IndelQDFilter")
    filteredIndels.filterExpression ++= List("\"QUAL<30.0\"","\"SB>-1.0\"","\"QD<2\"")
    filteredIndels.variantVCF = indels.out
    filteredIndels.out = swapExt("IndelCalls", indels.out, ".vcf",".filtered.vcf")

    // Call snps
    val snps = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    snps.analysisName = base+"_snps"
    snps.jobOutputFile = new File(".queue/logs/SNPCalling/UnifiedGenotyper.snps.out")
    snps.memoryLimit = Some(6)
    snps.downsample_to_coverage = Some(600)
    snps.input_file = bamFiles
    snps.rodBind :+= RodBind("dbsnp", "vcf", DBSNP)
    snps.out = new File("SnpCalls", base+".snps.vcf")

    snps.scatterCount = qscript.num_var_scatter_jobs
    snps.setupScatterFunction = {
      case scatter: ScatterFunction =>
        scatter.commandDirectory = new File("SnpCalls/ScatterGather")
        scatter.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Scatter.out")
    }
    snps.setupCloneFunction = {
      case (clone: CloneFunction, index: Int) =>
        clone.commandDirectory = new File("SnpCalls/ScatterGather/Scatter_%s".format(index))
        clone.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Scatter_%s.out".format(index))
    }
    snps.setupGatherFunction = {
      case (gather: GatherFunction, source: ArgumentSource) =>
        gather.commandDirectory = new File("SnpCalls/ScatterGather/Gather_%s".format(source.field.getName))
        gather.jobOutputFile = new File(".queue/logs/SNPCalling/ScatterGather/Gather_%s.out".format(source.field.getName))
    }

    // Filter snps at indels
    val filteredSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filteredSNPs.analysisName = base+"_filteredSNPs"
    filteredSNPs.jobOutputFile = new File(".queue/logs/SNPCalling/VariantFiltration.snps.out")
    filteredSNPs.clusterWindowSize = Some(10)
    filteredSNPs.clusterSize = Some(3)
    filteredSNPs.rodBind :+= RodBind("mask", "VCF", filteredIndels.out)
    filteredSNPs.variantVCF = snps.out
    filteredSNPs.out = swapExt("SnpCalls",snps.out,".vcf",".filtered.vcf")

    // Mako de Clusters
    val cr = new ContrastiveRecalibrator with CommandLineGATKArgs with ExpandedIntervals
    cr.rodBind :+= new RodBind("input","vcf",filteredSNPs.out)
    cr.rodBind :+= new RodBind("dbsnp","vcf",DBSNP,"known=true,training=false,truth=false,prior=8.0")
    cr.rodBind :+= new RodBind("hapmap","vcf", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"),"known=false,training=true,truth=true,prior=15.0")
    cr.rodBind :+= new RodBind("omni","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf"),"known=false,training=true,truth=true,prior=12.0")
    cr.allPoly = true
    cr.use_annotation ++= List("HaplotypeScore","SB","QD","HRun")
    cr.tranches_file = new File(base+".tranche")
    cr.recal_file = new File(base+".contrastive.recal.table")
    cr.tranche ++= List("99.9","99.5","99.25","98.0","97.75","97.65","97.5","97.3","97.2","97.1","98.0","97.5","97.0","96.75","96.5","96.0","95.5","95.0","94.75","94.5","94.25","94.0",
                        "93.75","93.5","93.25","93.0","92.75","92.5","92.25","92.0","91.0","90.0")
    cr.analysisName = base+"_ContrastiveRecalibrator"
    cr.memoryLimit = Some(32)
    cr.num_threads = Some(6)


    // Apply the Recalibration
    val ar = new ApplyRecalibration with CommandLineGATKArgs with ExpandedIntervals
    ar.rodBind :+= new RodBind("input","vcf",filteredSNPs.out)
    ar.tranches_file = cr.tranches_file
    ar.recal_file = cr.recal_file
    ar.ts_filter_level = Some(91.75)
    ar.out = new File(base+"_contrastive_recal.91.75.vcf")
    ar.memoryLimit = Some(6)

    // Variant eval the standard region
    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.analysisName = base+"_VariantEval"
    stdEval.jobOutputFile = new File(".queue/logs/Overall/VariantEval.std.out")
    stdEval.noST = true
    stdEval.noEV = true
    stdEval.evalModule ++= List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants","GenotypeConcordance")
    stdEval.stratificationModule ++= List("EvalRod", "CompRod", "Novelty","Sample")
    stdEval.rodBind :+= RodBind("dbsnp", "vcf",DBSNP)
    stdEval.rodBind :+= RodBind("evalContrastive", "VCF", ar.out)
    stdEval.rodBind :+= RodBind("evalHandFilter","VCF",handFiltered)
    stdEval.rodBind :+= RodBind("compHandFilter","VCF",handFiltered)
    stdEval.rodBind :+= RodBind("compAxiom","VCF",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/affymetrix_axiom/Affymetrix_Axiom_DB_2010_v4_b37.noOmni.noHM3.vcf"))
    stdEval.rodBind :+= RodBind("compOMNI","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf"))
    stdEval.out = swapExt(ar.out, ".vcf", ".eval")
    stdEval.num_threads = Some(6)

    // Variant eval the flanking region
    val flanksEval = new VariantEval with CommandLineGATKArgs
    flanksEval.analysisName = base+"_VariantEval"
    flanksEval.jobOutputFile = new File(".queue/logs/Overall/VariantEval.flanks.out")
    flanksEval.intervals = List(ei.outList)
    flanksEval.noST = true
    flanksEval.noEV = true
    flanksEval.evalModule ++= List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants","GenotypeConcordance")
    flanksEval.stratificationModule ++= List("EvalRod", "CompRod", "Novelty","Sample")
    flanksEval.rodBind :+= RodBind("dbsnp", "vcf",DBSNP)
    flanksEval.rodBind :+= RodBind("evalContrastive", "VCF", ar.out)
    flanksEval.rodBind :+= RodBind("evalHandFilter","VCF",handFiltered)
    flanksEval.rodBind :+= RodBind("compHandFilter","VCF",handFiltered)
    flanksEval.rodBind :+= RodBind("compAxiom","VCF",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/affymetrix_axiom/Affymetrix_Axiom_DB_2010_v4_b37.noOmni.noHM3.vcf"))
    flanksEval.rodBind :+= RodBind("compOMNI","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.vcf"))
    flanksEval.out = swapExt(ar.out, ".vcf", ".flanks.eval")
    flanksEval.num_threads = Some(6)

    // Make the bam list
    val listOfBams =  new File("Resources", base +".BamFiles.list")

    val writeBamList = new ListWriterFunction
    writeBamList.analysisName = base + "_BamList"
    writeBamList.jobOutputFile = new File(".queue/logs/Overall/WriteBamList.out")
    writeBamList.inputFiles = bamFiles
    writeBamList.listFile = listOfBams

    //add(indels, filteredIndels, snps, filteredSNPs, stdEval, writeBamList,cr,ar)
    add(ar,stdEval)
    
    if (qscript.expandIntervals > 0) {
      add(flanksEval)
    }

  }
}
