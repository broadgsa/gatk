
import java.io.PrintWriter
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class fullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="path to trigger track (for UnifiedGenotyper)", shortName="trigger", required=false)
  var trigger: File = _

  @Input(doc="path to refseqTable (for GenomicAnnotator)", shortName="refseqTable")
  var refseqTable: File = _

  @Input(doc="path to Picard FixMateInformation.jar.  See http://picard.sourceforge.net/ .", required=false)
  var picardFixMatesJar: File = new java.io.File("/seq/software/picard/current/bin/FixMateInformation.jar")

  @Input(doc="path to GATK jar")
  var gatkJar: File = _

  @Input(doc="per-sample downsampling level",shortName="dcov",required=false)
  var downsampling_coverage = 600

  @Input(doc="level of parallelism for IndelRealigner.  By default is set to 1.", shortName="cleanerScatter", required=false)
  var num_cleaner_scatter_jobs = 1

  @Input(doc="level of parallelism for UnifiedGenotyper.   By default is set to 20.", shortName="snpScatter", required=false)
  var num_snp_scatter_jobs = 20

  //@Input(doc="level of parallelism for IndelGenotyperV2", shortName="indelScatter", required=false)
  //var num_indel_scatter_jobs = 5

  @Input(doc="Skip indel-cleaning for BAM files (for testing only)", shortName="skipCleaning", required=false)
  var skip_cleaning = false

  @Input(doc="ADPR script", shortName ="tearScript", required=true)
  var tearScript: File = _

  //@Input(doc="Sequencing maching name (for use by adpr)")
  //var machine: String = _

  //@Input(doc="Sequencing experiement type (for use by adpr)--Whole_Exome, Whole_Genome, or Hybrid_Selection")
  //var protocol: String = _

  @Argument(doc="Job queue for large memory jobs (>4 to 16GB)", shortName="bigMemQueue", required=false)
  var big_mem_queue: String = _

  @Argument(doc="Job queue for short run jobs (<1hr)", shortName="shortJobQueue", required=false)
  var short_job_queue: String = _

  private var pipeline: Pipeline = _

  private var dbsnpType: String = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals = List(qscript.pipeline.getProject.getIntervalList)
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.memoryLimit = Some(4)
  }


  // ------------ SETUP THE PIPELINE ----------- //


  def script = {
    pipeline = YamlUtils.load(classOf[Pipeline], qscript.yamlFile)
    //var dbsnpType: String = _ //figure out how to get it so this is recognized and the whole thing has access to it.
      if (qscript.pipeline.getProject.getDbsnpFile.toString.contains("rod")){
        dbsnpType = "dbSNP"
      } else  {
        dbsnpType = "VCF"
      }


    val projectBase: String = qscript.pipeline.getProject.getName
    if (qscript.skip_cleaning) {
      //endToEnd(projectBase + ".uncleaned", "recalibrated", adprRscript, seq, expKind)

      endToEnd(projectBase + ".uncleaned", "recalibrated")
    } else {

      // there are commands that use all the bam files
      val recalibratedSamples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains("recalibrated"))
      //val adprRScript = qscript.adprScript
      //val seq = qscript.machine
      //val expKind = qscript.protocol

      // get max num contigs for indel cleaning parallelism, plus 1 for -L unmapped
      val numContigs = IntervalUtils.distinctContigs(qscript.pipeline.getProject.getReferenceFile).size + 1

      for ( sample <- recalibratedSamples ) {
        val sampleId = sample.getId
        // put unclean bams in unclean genotypers in advance, create the extension files
        val bam = sample.getBamFiles.get("recalibrated")
        if (!sample.getBamFiles.contains("cleaned")) {
          sample.getBamFiles.put("cleaned", swapExt("CleanedBams", bam,"bam","cleaned.bam"))
        }

        val cleaned_bam = sample.getBamFiles.get("cleaned")
        val indel_targets = swapExt("CleanedBams/IntermediateFiles/"+sampleId, bam,"bam","realigner_targets.interval_list")

        // create the cleaning commands
        val targetCreator = new RealignerTargetCreator with CommandLineGATKArgs
        targetCreator.jobOutputFile = new File(".queue/logs/Cleaning/%s/RealignerTargetCreator.out".format(sampleId))
        targetCreator.jobName = "CreateTargets_"+sampleId
        targetCreator.analysisName = "CreateTargets_"+sampleId
        targetCreator.input_file :+= bam
        targetCreator.out = indel_targets
        targetCreator.memoryLimit = Some(2)
        targetCreator.isIntermediate = true

        val realigner = new IndelRealigner with CommandLineGATKArgs
        realigner.jobOutputFile = new File(".queue/logs/Cleaning/%s/IndelRealigner.out".format(sampleId))
        realigner.analysisName = "RealignBam_"+sampleId
        realigner.input_file = targetCreator.input_file
        realigner.targetIntervals = targetCreator.out
        realigner.intervals = Nil
        realigner.intervalsString = Nil
        realigner.scatterCount = num_cleaner_scatter_jobs min numContigs
        realigner.rodBind :+= RodBind("dbsnp", dbsnpType, qscript.pipeline.getProject.getDbsnpFile)
        realigner.rodBind :+= RodBind("indels", "VCF", swapExt(realigner.reference_sequence.getParentFile, realigner.reference_sequence, "fasta", "1kg_pilot_indels.vcf"))

        // if scatter count is > 1, do standard scatter gather, if not, explicitly set up fix mates
        if (realigner.scatterCount > 1) {
          realigner.out = cleaned_bam
          // While gathering run fix mates.
          realigner.setupScatterFunction = {
            case scatter: ScatterFunction =>
              scatter.commandDirectory = new File("CleanedBams/IntermediateFiles/%s/ScatterGather".format(sampleId))
              scatter.jobOutputFile = new File(".queue/logs/Cleaning/%s/Scatter.out".format(sampleId))
          }
          realigner.setupCloneFunction = {
            case (clone: CloneFunction, index: Int) =>
              clone.commandDirectory = new File("CleanedBams/IntermediateFiles/%s/ScatterGather/Scatter_%s".format(sampleId, index))
              clone.jobOutputFile = new File(".queue/logs/Cleaning/%s/Scatter_%s.out".format(sampleId, index))
          }
          realigner.setupGatherFunction = {
            case (gather: BamGatherFunction, source: ArgumentSource) =>
              gather.commandDirectory = new File("CleanedBams/IntermediateFiles/%s/ScatterGather/Gather_%s".format(sampleId, source.field.getName))
              gather.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
              gather.memoryLimit = Some(6)
              gather.jarFile = qscript.picardFixMatesJar
              // Don't pass this AS=true to fix mates!
              gather.assumeSorted = None
            case (gather: GatherFunction, source: ArgumentSource) =>
              gather.commandDirectory = new File("CleanedBams/IntermediateFiles/%s/ScatterGather/Gather_%s".format(sampleId, source.field.getName))
              gather.jobOutputFile = new File(".queue/logs/Cleaning/%s/Gather_%s.out".format(sampleId, source.field.getName))
          }

          add(targetCreator,realigner)
        } else {
          realigner.out = swapExt("CleanedBams/IntermediateFiles/"+sampleId,bam,"bam","unfixed.cleaned.bam")
          realigner.isIntermediate = true

          // Explicitly run fix mates if the function won't be scattered.
          val fixMates = new PicardBamJarFunction {
            // Declare inputs/outputs for dependency tracking.
            @Input(doc="unfixed bam") var unfixed: File = _
            @Output(doc="fixed bam") var fixed: File = _
            def inputBams = List(unfixed)
            def outputBam = fixed
          }

          fixMates.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
          fixMates.memoryLimit = Some(6)
          fixMates.jarFile = qscript.picardFixMatesJar
          fixMates.unfixed = realigner.out
          fixMates.fixed = cleaned_bam
          fixMates.analysisName = "FixMates_"+sampleId

          // Add the fix mates explicitly
          add(targetCreator,realigner,fixMates)
        }

        var samtoolsindex = new SamtoolsIndexFunction
        samtoolsindex.jobOutputFile = new File(".queue/logs/Cleaning/%s/SamtoolsIndex.out".format(sampleId))
        samtoolsindex.bamFile = cleaned_bam
        samtoolsindex.analysisName = "index_cleaned_"+sampleId
        samtoolsindex.jobQueue = qscript.short_job_queue
        add(samtoolsindex)
      }

      //endToEnd(projectBase + ".cleaned", "cleaned", adprRscript, seq, expKind)
      endToEnd(projectBase + ".cleaned", "cleaned")
    }
  }

  def endToEnd(base: String, bamType: String) = {

    val samples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains(bamType)).toList
    val bamFiles = samples.map(_.getBamFiles.get(bamType))

    // step through the un-indel-cleaned graph:
    // 1a. call snps and indels
    val snps = new UnifiedGenotyper with CommandLineGATKArgs
    snps.jobOutputFile = new File(".queue/logs/SNPCalling/UnifiedGenotyper.out")
    snps.analysisName = base+"_SNP_calls"
    snps.input_file = bamFiles
    //snps.annotation ++= List("AlleleBalance")
    snps.input_file = bamFiles
    snps.group :+= "Standard"
    snps.out = new File("SnpCalls", base+".vcf")
    //snps.standard_min_confidence_threshold_for_emitting = Some(10)
    //snps.min_mapping_quality_score = Some(20)
    //snps.min_base_quality_score = Some(20)
    snps.downsample_to_coverage = Some(qscript.downsampling_coverage)
    //snps.annotation :+= "QualByDepthV2"
    snps.rodBind :+= RodBind("dbsnp", dbsnpType, qscript.pipeline.getProject.getDbsnpFile)
    snps.memoryLimit = Some(6)
    //if (qscript.trigger != null) {
    //  snps.trigger_min_confidence_threshold_for_calling = Some(30)
    //  snps.rodBind :+= RodBind("trigger", "VCF", qscript.trigger)
    //  // TODO: triggers need to get a name for comp-ing them if not dbSNP?
    //  snps.rodBind :+= RodBind( "compTrigger", "VCF", qscript.trigger )
    //}

    // todo -- add generalize comp inputs
    //if ( qscript.comp1KGCEU != null ) {
    //  snps.rodBind :+= RodBind( "comp1KG_CEU", "VCF", qscript.comp1KGCEU )
    //}

    snps.scatterCount = qscript.num_snp_scatter_jobs
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

    // indel genotyper does one sample at a time
    var indelCallFiles = List.empty[RodBind]
    var indelGenotypers = List.empty[IndelGenotyperV2 with CommandLineGATKArgs]
    var loopNo = 0
    var priority = ""
    for ( sample <- samples ) {
      val sampleId = sample.getId
      val bam = sample.getBamFiles.get(bamType)

      var indel = new IndelGenotyperV2 with CommandLineGATKArgs
      indel.jobOutputFile = new File(".queue/logs/IndelCalling/%s/IndelGenotyperV2.out".format(sampleId))
      indel.window_size = Some(350)
      indel.analysisName = "IndelGenotyper_"+sampleId
      indel.input_file :+= bam
      indel.out = swapExt("IndelCalls/IntermediateFiles/" + sampleId, bam,".bam",".indels.vcf")
      indel.downsample_to_coverage = Some(qscript.downsampling_coverage)
      indelCallFiles :+= RodBind("v"+loopNo.toString, "VCF", indel.out)
      //indel.scatterCount = qscript.num_indel_scatter_jobs

      indelGenotypers :+= indel

      if ( loopNo == 0 ) {
        priority = "v0"
      } else {
        priority += ",v"+loopNo.toString
      }
      loopNo += 1
    }
    val mergeIndels = new CombineVariants with CommandLineGATKArgs
    mergeIndels.jobOutputFile = new File(".queue/logs/IndelCalling/CombineVariants.out")
    mergeIndels.out = new TaggedFile("IndelCalls/" + qscript.pipeline.getProject.getName+".indels.vcf","vcf")
    mergeIndels.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNSORTED)
    mergeIndels.priority = priority
    mergeIndels.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    mergeIndels.rodBind = indelCallFiles
    mergeIndels.analysisName = base+"_MergeIndels"
    mergeIndels.memoryLimit = Some(4)
    mergeIndels.jobQueue = qscript.big_mem_queue

    // 1b. genomically annotate SNPs -- no longer slow
    val annotated = new GenomicAnnotator with CommandLineGATKArgs
    annotated.jobOutputFile = new File(".queue/logs/SNPCalling/GenomicAnnotator.out")
    annotated.rodBind :+= RodBind("variant", "VCF", snps.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.refseqTable)
    //annotated.rodBind :+= RodBind("dbsnp", "AnnotatorInputTable", qscript.dbsnpTable)
    annotated.out = swapExt("SnpCalls",snps.out,".vcf",".annotated.vcf")
    //annotated.select :+= "dbsnp.name,dbsnp.refUCSC,dbsnp.strand,dbsnp.observed,dbsnp.avHet"
    annotated.rodToIntervalTrackName = "variant"
    annotated.analysisName = base+"_GenomicAnnotator"

    // 2.a filter on cluster and near indels
    val masker = new VariantFiltration with CommandLineGATKArgs
    masker.jobOutputFile = new File(".queue/logs/SNPCalling/Masker.out")
    masker.variantVCF = annotated.out
    masker.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    masker.maskName = "NearIndel"
    masker.clusterWindowSize = Some(10)
    masker.clusterSize = Some(3)
    masker.out = swapExt("SnpCalls",annotated.out,".vcf",".indel.masked.vcf")
    masker.analysisName = base+"_Cluster_and_Indel_filter"

    // 2.b hand filter with standard filter
    val handFilter = new VariantFiltration with CommandLineGATKArgs
    handFilter.jobOutputFile = new File(".queue/logs/SNPCalling/HandFilter.out")
    handFilter.variantVCF = masker.out
    handFilter.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    //handFilter.filterName ++= List("StrandBias","AlleleBalance","QualByDepth","HomopolymerRun")
    //handFilter.filterExpression ++= List("\"SB>=0.10\"","\"AB>=0.75\"","\"QD<5.0\"","\"HRun>=4\"")
    handFilter.filterName ++= List("StrandBias","QualByDepth","HomopolymerRun")
    handFilter.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    handFilter.out = swapExt("SnpCalls",annotated.out,".vcf",".handfiltered.vcf")
    handFilter.analysisName = base+"_HandFilter"

    // 3.i generate gaussian clusters on the masked vcf
    // todo -- args for annotations?
    // todo -- args for resources (properties file)
   // val clusters = new GenerateVariantClusters with CommandLineGATKArgs
  //  clusters.jobOutputFile = new File(".queue/logs/SNPCalling/Clusters.out")
  //  clusters.rodBind :+= RodBind("input", "VCF", masker.out)
  //  clusters.rodBind :+= RodBind("dbsnp", "ROD", qscript.pipeline.getProject.getDbsnpFile)
  //  val clusters_clusterFile = swapExt("SnpCalls/IntermediateFiles",snps.out,".vcf",".cluster")
  //  clusters.clusterFile = clusters_clusterFile
  //  clusters.memoryLimit = Some(4)
  //  clusters.jobQueue = qscript.big_mem_queue

   // clusters.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
   // clusters.analysisName = base+"_Cluster"

    // 3.ii apply gaussian clusters to the masked vcf
   // val recalibrate = new VariantRecalibrator with CommandLineGATKArgs
   // recalibrate.jobOutputFile = new File(".queue/logs/SNPCalling/Recalibrator.out")
   // recalibrate.clusterFile = clusters.clusterFile
   // recalibrate.DBSNP = qscript.pipeline.getProject.getDbsnpFile
  //  recalibrate.rodBind :+= RodBind("input", "VCF", masker.out)
  //  recalibrate.out = swapExt("SnpCalls",masker.out,".vcf",".recalibrated.vcf")
//    recalibrate.target_titv = qscript.target_titv
   // recalibrate.tranches_file = swapExt("SnpCalls/IntermediateFiles", masker.out,".vcf",".recalibrate.tranches")
   // recalibrate.analysisName = base+"_VariantRecalibrator"

    // 3.iii apply variant cuts to the clusters
   // val cut = new ApplyVariantCuts with CommandLineGATKArgs
   // cut.jobOutputFile = new File(".queue/logs/SNPCalling/VariantCuts.out")
   // cut.rodBind :+= RodBind("input", "VCF", recalibrate.out)
   // cut.out = swapExt("SnpCalls",recalibrate.out,".vcf",".tranched.vcf")
    //cut.tranches_file = recalibrate.tranches_file
    // todo -- fdr inputs, etc
//    cut.fdr_filter_level = Some(1)
  //  cut.analysisName = base+"_ApplyVariantCuts"

    // 4. Variant eval the cut and the hand-filtered vcf files
    val eval = new VariantEval with CommandLineGATKArgs
    eval.jobOutputFile = new File(".queue/logs/SNPCalling/VariantEval.out")
   // eval.rodBind :+= RodBind("evalOptimized", "VCF", cut.out)
    eval.rodBind :+= RodBind("eval", "VCF", handFilter.out)
    eval.evalModule ++= List("FunctionalClassBySample", "SimpleMetricsBySample", "CountFunctionalClasses", "CompOverlap", "CountVariants", "TiTvVariantEvaluator")
    eval.reportLocation = new File("SnpCalls", base+".eval")
    eval.reportType = Option(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.R)
    eval.analysisName = base+"_VariantEval"
    if(dbsnpType=="VCF"){
     eval.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")

    } else{
      eval.rodBind :+= RodBind("dbsnp", dbsnpType, qscript.pipeline.getProject.getDbsnpFile)
    }
    add(snps)

    // 5. Make the bam list
    val listOfBams =  new File(base +".BamFiles.list")

    class BamListWriter extends InProcessFunction {
      @Input(doc="bamFiles") var bamFiles: List[File] = Nil
      @Output(doc="bamList") var bamList: File = _

      def run {
        val writer = new PrintWriter(bamList)
        for (bamFile <- bamFiles)
          writer.println(bamFile.toString)
        writer.close()
      }
    }

    val writeBamList = new BamListWriter
    writeBamList.bamFiles = bamFiles
    writeBamList.bamList = listOfBams
    writeBamList.analysisName = base + "_BamList"
    writeBamList.jobOutputFile = new File(".queue/logs/SNPCalling/bamlist.out")

    // 6. Run the ADPR and make pretty stuff

    class rCommand extends CommandLineFunction{
      @Argument(doc="R script")
      var  script: File = _
      @Argument(doc="list of bams")
        var bamlist: File =_
      @Input(doc="Eval files root")
      var evalroot: File =_
      @Output(doc="plot pdf loc")
      var plots: File =_
      @Output(doc="tearsheet loc")
      var tearsheet: File =_
      def commandLine = "Rscript %s -bamlist %s -evalroot %s -tearout %s -plotout %s"
        .format(script, bamlist, evalroot, tearsheet, plots)
    }

    val adpr = new rCommand
     adpr.bamlist = listOfBams
     adpr.script = tearScript
     adpr.evalroot = eval.reportLocation
     adpr.jobOutputFile = new File(".queue/logs/SNPCalling/adpr.out")
     adpr.tearsheet = new File("SnpCalls", base + ".tearsheet.pdf")
     adpr.plots = new File("SnpCalls", base + ".plots.pdf")
     adpr.analysisName = base + "_ADPR"


    for ( igv2 <- indelGenotypers ) {
      add(igv2)
    }

    add(mergeIndels,annotated,masker,handFilter,eval,writeBamList,adpr)
  }
}
