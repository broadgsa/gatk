import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.broad.PicardPipeline

class FullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="path to GATK jar", shortName="G")
  var gatkJar: File = _

  @Input(doc="level of parallelism for IndelRealigner.  By default is set to 1.", shortName="cleanerScatter", required=false)
  var num_cleaner_scatter_jobs = 1

  @Input(doc="level of parallelism for UnifiedGenotyper (both for SNPs and indels).  By default is set to 20.", shortName="varScatter", required=false)
  var num_var_scatter_jobs = 20

  @Argument(doc="expand each target in input intervals by the specified number of bases (50 bases by default)", shortName="expand", required=false)
  var expandIntervals = 50

  @Input(doc="Skip indel-cleaning for BAM files (for testing only)", shortName="skipCleaning", required=false)
  var skip_cleaning = false

  @Input(doc="ADPR script", shortName ="tearScript", required=false)
  var tearScript: File = _

  private var pipeline: Pipeline = _

  private final val picardFixMatesClass = "net.sf.picard.sam.FixMateInformation"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals = List(qscript.pipeline.getProject.getIntervalList)
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.memoryLimit = 4
  }


  // ------------ SETUP THE PIPELINE ----------- //


  def script() {
    pipeline = PicardPipeline.parse(qscript.yamlFile)

    val projectBase: String = qscript.pipeline.getProject.getName

    if (qscript.skip_cleaning) {
      endToEnd(projectBase + ".uncleaned", "recalibrated")
    } else {
      val recalibratedSamples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains("recalibrated"))

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
        targetCreator.memoryLimit = 2
        targetCreator.isIntermediate = true

        val realigner = new IndelRealigner with CommandLineGATKArgs
        realigner.jobOutputFile = new File(".queue/logs/Cleaning/%s/IndelRealigner.out".format(sampleId))
        realigner.analysisName = "RealignBam_"+sampleId
        realigner.input_file = targetCreator.input_file
        realigner.targetIntervals = targetCreator.out
        realigner.intervals = Nil
        realigner.intervalsString = Nil
        realigner.scatterCount = num_cleaner_scatter_jobs
        realigner.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
        realigner.rodBind :+= RodBind("indels", "VCF", swapExt(realigner.reference_sequence.getParentFile, realigner.reference_sequence, "fasta", "1kg_pilot_indels.vcf"))

        // if scatter count is > 1, do standard scatter gather, if not, explicitly set up fix mates
        if (realigner.scatterCount > 1) {
          realigner.out = cleaned_bam
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
              gather.memoryLimit = 6
              gather.javaMainClass = picardFixMatesClass
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
          val fixMates = new PicardBamFunction {
            @Input(doc="unfixed bam") var unfixed: File = _
            @Output(doc="fixed bam") var fixed: File = _
            def inputBams = List(unfixed)
            def outputBam = fixed
          }

          fixMates.jobOutputFile = new File(".queue/logs/Cleaning/%s/FixMates.out".format(sampleId))
          fixMates.memoryLimit = 6
          fixMates.javaMainClass = picardFixMatesClass
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
        //samtoolsindex.jobQueue = qscript.short_job_queue
        add(samtoolsindex)
      }

      endToEnd(projectBase + ".cleaned", "cleaned")
    }
  }

  def endToEnd(base: String, bamType: String) = {
    val samples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains(bamType)).toList
    val bamFiles = samples.map(_.getBamFiles.get(bamType))

    val ei : ExpandIntervals = new ExpandIntervals(qscript.pipeline.getProject.getIntervalList, 1, qscript.expandIntervals, new File("Resources", base + ".flanks.interval_list"), qscript.pipeline.getProject.getReferenceFile, "INTERVALS", "INTERVALS")
    ei.jobOutputFile = new File(".queue/logs/Overall/ExpandIntervals.out")

    if (qscript.expandIntervals > 0) {
      add(ei)
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
    indels.memoryLimit = 6
    indels.downsample_to_coverage = 600
    indels.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    indels.input_file = bamFiles
    indels.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
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
    snps.memoryLimit = 6
    snps.downsample_to_coverage = 600
    snps.input_file = bamFiles
    snps.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    snps.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getGenotypeDbsnpType, qscript.pipeline.getProject.getGenotypeDbsnp)
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

    // Filter snps
    val filteredSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filteredSNPs.analysisName = base+"_filteredSNPs"
    filteredSNPs.jobOutputFile = new File(".queue/logs/SNPCalling/VariantFiltration.snps.out")
    filteredSNPs.filterName ++= List("SNPSBFilter","SNPQDFilter","SNPHRunFilter")
    filteredSNPs.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    filteredSNPs.clusterWindowSize = 10
    filteredSNPs.clusterSize = 3
    filteredSNPs.rodBind :+= RodBind("mask", "VCF", filteredIndels.out)
    filteredSNPs.variantVCF = snps.out
    filteredSNPs.out = swapExt("SnpCalls",snps.out,".vcf",".filtered.vcf")

    // Combine indels and snps into one VCF
    val combineAll = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineAll.analysisName = base + "_combineAll"
    combineAll.jobOutputFile = new File(".queue/logs/Combined/CombineVariants.out")
    combineAll.variantMergeOptions = org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION
    combineAll.rod_priority_list = "Indels,SNPs"
    combineAll.rodBind :+= RodBind("Indels", "VCF", filteredIndels.out)
    combineAll.rodBind :+= RodBind("SNPs", "VCF", filteredSNPs.out)
    combineAll.out = new File("CombinedCalls", base + ".allVariants.filtered.vcf")

    // Annotate variants
    val annotated = new GenomicAnnotator with CommandLineGATKArgs with ExpandedIntervals
    annotated.analysisName = base+"_annotated"
    annotated.jobOutputFile = new File(".queue/logs/Combined/GenomicAnnotator.out")
    annotated.rodToIntervalTrackName = "variant"
    annotated.rodBind :+= RodBind("variant", "VCF", combineAll.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.pipeline.getProject.getRefseqTable)
    annotated.out = new File(base + ".snps_and_indels.filtered.annotated.vcf")

    // Variant eval the standard region
    val stdEval = new VariantEval with CommandLineGATKArgs
    stdEval.analysisName = base+"_VariantEval"
    stdEval.jobOutputFile = new File(".queue/logs/Overall/VariantEval.std.out")
    stdEval.noST = true
    stdEval.noEV = true
    stdEval.evalModule ++= List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    stdEval.stratificationModule ++= List("EvalRod", "CompRod", "Novelty")
    stdEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
    stdEval.rodBind :+= RodBind("eval", "VCF", annotated.out)
    stdEval.out = swapExt(annotated.out, ".vcf", ".eval")

    // Variant eval the flanking region
    val flanksEval = new VariantEval with CommandLineGATKArgs
    flanksEval.analysisName = base+"_VariantEval"
    flanksEval.jobOutputFile = new File(".queue/logs/Overall/VariantEval.flanks.out")
    flanksEval.intervals = List(ei.outList)
    flanksEval.noST = true
    flanksEval.noEV = true
    flanksEval.evalModule ++= List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    flanksEval.stratificationModule ++= List("EvalRod", "CompRod", "Novelty")
    flanksEval.rodBind :+= RodBind("dbsnp", qscript.pipeline.getProject.getEvalDbsnpType, qscript.pipeline.getProject.getEvalDbsnp)
    flanksEval.rodBind :+= RodBind("eval", "VCF", annotated.out)
    flanksEval.out = swapExt(annotated.out, ".vcf", ".flanks.eval")

    // Make the bam list
    val listOfBams =  new File("Resources", base +".BamFiles.list")

    val writeBamList = new ListWriterFunction
    writeBamList.analysisName = base + "_BamList"
    writeBamList.jobOutputFile = new File(".queue/logs/Overall/WriteBamList.out")
    writeBamList.inputFiles = bamFiles
    writeBamList.listFile = listOfBams

    add(indels, filteredIndels, snps, filteredSNPs, combineAll, annotated, stdEval, writeBamList)
    
    if (qscript.expandIntervals > 0) {
      add(flanksEval)
    }

    // Run the ADPR and make pretty stuff
    if (qscript.tearScript != null) {
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
     adpr.analysisName = base + "_ADPR"
     adpr.bamlist = listOfBams
     adpr.yaml = qscript.yamlFile.getAbsoluteFile
     adpr.script = qscript.tearScript
     adpr.evalroot = stdEval.out
     adpr.jobOutputFile = new File(".queue/logs/Overall/ADPR.out")
     adpr.tearsheet = new File("VariantCalls", base + ".tearsheet.pdf")

     add(adpr)
    }
  }
}
