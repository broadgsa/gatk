
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils

class fullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="path to trigger track (for UnifiedGenotyper)", shortName="trigger", required=false)
  var trigger: File = _

  // TODO: Fix command lines that pass -refseqTable
  @Input(doc="path to refseqTable (for GenomicAnnotator) if not present in the YAML", shortName="refseqTable", required=false)
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

  @Input(doc="ADPR script", shortName ="tearScript", required=false)
  var tearScript: File = _

  //@Input(doc="Sequencing maching name (for use by adpr)")
  //var machine: String = _

  //@Input(doc="Sequencing experiement type (for use by adpr)--Whole_Exome, Whole_Genome, or Hybrid_Selection")
  //var protocol: String = _

  // TODO: Fix command lines that pass -bigMemQueue
  @Argument(doc="Unused", shortName="bigMemQueue", required=false)
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
    // TODO: Fix command lines that pass -refseqTable
    if (qscript.refseqTable != null)
      qscript.pipeline.getProject.setRefseqTable(qscript.refseqTable)
    if (qscript.skip_cleaning) {
      //endToEnd(projectBase + ".uncleaned", "recalibrated", adprRscript, seq, expKind)

      endToEnd(projectBase + ".uncleaned", "recalibrated")
    } else {

      // there are commands that use all the bam files
      val recalibratedSamples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains("recalibrated"))
      //val adprRScript = qscript.adprScript
      //val seq = qscript.machine
      //val expKind = qscript.protocol

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
        realigner.scatterCount = num_cleaner_scatter_jobs
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
    snps.input_file = bamFiles
    snps.group :+= "Standard"
    snps.out = new File("SnpCalls", base+".vcf")
    snps.downsample_to_coverage = Some(qscript.downsampling_coverage)
    snps.rodBind :+= RodBind("dbsnp", dbsnpType, qscript.pipeline.getProject.getDbsnpFile)
    snps.memoryLimit = Some(6)

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

    val indels = new UnifiedGenotyper with CommandLineGATKArgs
    indels.jobOutputFile = new File(".queue/logs/IndelCalling/UnifiedGenotyper.out")
    indels.analysisName = base+"_Indel_calls"
    indels.input_file = bamFiles
    indels.input_file = bamFiles
    indels.group :+= "Standard"
    indels.out = new File("IndelCalls", base+".vcf")
    indels.downsample_to_coverage = Some(qscript.downsampling_coverage)
    indels.rodBind :+= RodBind("dbsnp", dbsnpType, qscript.pipeline.getProject.getDbsnpFile)
    indels.memoryLimit = Some(6)
    indels.genotype_likelihoods_model = Option(org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.DINDEL)



    indels.scatterCount = qscript.num_snp_scatter_jobs
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


    // 1b. genomically annotate SNPs -- no longer slow
    val annotated = new GenomicAnnotator with CommandLineGATKArgs
    annotated.jobOutputFile = new File(".queue/logs/SNPCalling/GenomicAnnotator.out")
    annotated.rodBind :+= RodBind("variant", "VCF", snps.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.pipeline.getProject.getRefseqTable)
    annotated.out = swapExt("SnpCalls",snps.out,".vcf",".annotated.vcf")
    annotated.rodToIntervalTrackName = "variant"
    annotated.analysisName = base+"_GenomicAnnotator"

    // 2.a filter on cluster and near indels
    val masker = new VariantFiltration with CommandLineGATKArgs
    masker.jobOutputFile = new File(".queue/logs/SNPCalling/Masker.out")
    masker.variantVCF = annotated.out
    masker.rodBind :+= RodBind("mask", "VCF", indels.out)
    masker.maskName = "NearIndel"
    masker.clusterWindowSize = Some(10)
    masker.clusterSize = Some(3)
    masker.out = swapExt("SnpCalls",annotated.out,".vcf",".indel.masked.vcf")
    masker.analysisName = base+"_Cluster_and_Indel_filter"

    // 2.b hand filter with standard filter
    val handFilter = new VariantFiltration with CommandLineGATKArgs
    handFilter.jobOutputFile = new File(".queue/logs/SNPCalling/HandFilter.out")
    handFilter.variantVCF = masker.out
    handFilter.rodBind :+= RodBind("mask", "VCF", indels.out)
    //handFilter.filterName ++= List("StrandBias","AlleleBalance","QualByDepth","HomopolymerRun")
    //handFilter.filterExpression ++= List("\"SB>=0.10\"","\"AB>=0.75\"","\"QD<5.0\"","\"HRun>=4\"")
    handFilter.filterName ++= List("StrandBias","QualByDepth","HomopolymerRun")
    handFilter.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    handFilter.out = swapExt("SnpCalls",annotated.out,".vcf",".handfiltered.vcf")
    handFilter.analysisName = base+"_HandFilter"


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

    // 5. Make the bam list
    val listOfBams =  new File(base +".BamFiles.list")

    val writeBamList = new ListWriterFunction
    writeBamList.inputFiles = bamFiles
    writeBamList.listFile = listOfBams
    writeBamList.analysisName = base + "_BamList"
    writeBamList.jobOutputFile = new File(".queue/logs/SNPCalling/bamlist.out")

    // 6. Run the ADPR and make pretty stuff

    add(snps, indels, annotated,masker,handFilter,eval,writeBamList)
    if (qscript.tearScript != null){
      class rCommand extends CommandLineFunction{
        @Input(doc="R script")
        var script: File = _
        @Input(doc="pipeline yaml")
        var yaml: File = _
        @Input(doc="list of bams")
         var bamlist: File =_
        @Input(doc="Eval files root")
        var evalroot: File =_
        @Output(doc="tearsheet loc")
        var tearsheet: File =_
        def commandLine = "Rscript %s -yaml %s -bamlist %s -evalroot %s -tearout %s".format(script, yaml, bamlist, evalroot, tearsheet)
      }

     val adpr = new rCommand
     adpr.bamlist = listOfBams
     adpr.yaml = qscript.yamlFile.getAbsoluteFile
     adpr.script = qscript.tearScript
     adpr.evalroot = eval.reportLocation
     adpr.jobOutputFile = new File(".queue/logs/SNPCalling/adpr.out")
     adpr.tearsheet = new File("SnpCalls", base + ".tearsheet.pdf")
     adpr.analysisName = base + "_ADPR"


     add(adpr)
    }




  }
}
