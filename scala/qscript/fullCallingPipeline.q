import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class fullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="list of contigs in the reference over which indel-cleaning jobs should be scattered (ugly)", shortName="contigIntervals")
  var contigIntervals: File = _

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

  @Input(doc="target Ti/Tv ratio for recalibration", shortName="titv", required=true)
  var target_titv: Float = _

  @Input(doc="per-sample downsampling level",shortName="dcov",required=false)
  var downsampling_coverage = 300

  @Input(doc="level of parallelism for UnifiedGenotyper", shortName="snpScatter", required=false)
  var num_snp_scatter_jobs = 20

  @Input(doc="level of parallelism for IndelGenotyperV2", shortName="indelScatter", required=false)
  var num_indel_scatter_jobs = 5

  @Input(doc="Skip indel-cleaning for BAM files (for testing only)", shortName="skipCleaning", required=false)
  var skip_cleaning = false

  //@Input(doc="ADPR script")
  //var adprScript: File = _

  //@Input(doc="Sequencing maching name (for use by adpr)")
  //var machine: String = _

  //@Input(doc="Sequencing experiement type (for use by adpr)--Whole_Exome, Whole_Genome, or Hybrid_Selection")
  //var protocol: String = _

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals = qscript.pipeline.getProject.getIntervalList
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.memoryLimit = Some(6)
  }


  // ------------ SETUP THE PIPELINE ----------- //


  def script = {
    pipeline = YamlUtils.load(classOf[Pipeline], qscript.yamlFile)

    val projectBase: String = qscript.pipeline.getProject.getName
    val cleanedBase: String = projectBase + ".cleaned"
    val uncleanedBase: String = projectBase + ".uncleaned"

    // there are commands that use all the bam files
    val recalibratedSamples = qscript.pipeline.getSamples.filter(_.getBamFiles.contains("recalibrated"))
    //val adprRScript = qscript.adprScript
    //val seq = qscript.machine
    //val expKind = qscript.protocol

    // count number of contigs (needed for indel cleaning parallelism)
    var contigCount = 0
    for ( line <- scala.io.Source.fromFile(qscript.contigIntervals).getLines ) {
      contigCount += 1
    }

    for ( sample <- recalibratedSamples ) {
      // put unclean bams in unclean genotypers in advance, create the extension files
      val bam = sample.getBamFiles.get("recalibrated")
      if (!sample.getBamFiles.contains("cleaned")) {
        sample.getBamFiles.put("cleaned", swapExt(bam,"bam","cleaned.bam"))
      }

      val cleaned_bam = sample.getBamFiles.get("cleaned")
      val indel_targets = swapExt(bam,"bam","realigner_targets.interval_list")

      // create the cleaning commands
      val targetCreator = new RealignerTargetCreator with CommandLineGATKArgs
      targetCreator.analysisName = "CreateTargets_"+bam.getName
      targetCreator.input_file :+= bam
      targetCreator.out = indel_targets

      val realigner = new IndelRealigner with CommandLineGATKArgs
      realigner.analysisName = "RealignBam_"+bam.getName
      realigner.input_file = targetCreator.input_file
      realigner.intervals = qscript.contigIntervals
      realigner.targetIntervals = new java.io.File(targetCreator.out.getAbsolutePath)
      realigner.scatterCount = contigCount

      // may need to explicitly run fix mates
      var fixMates = new PicardBamJarFunction {
          // Declare inputs/outputs for dependency tracking.
          @Input(doc="unfixed bam") var unfixed: File = _
          @Output(doc="fixed bam") var fixed: File = _
          def inputBams = List(unfixed)
          def outputBam = fixed
        }

      // realigner.out = cleaned_bam
      // realigner.scatterClass = classOf[ContigScatterFunction]
      // realigner.setupGatherFunction = { case (f: BamGatherFunction, _) => f.jarFile = qscript.picardFixMatesJar }
      // realigner.jobQueue = "week"

      // if scatter count is > 1, do standard scatter gather, if not, explicitly set up fix mates
      if (realigner.scatterCount > 1) {
        realigner.out = cleaned_bam
        // While gathering run fix mates.
        realigner.scatterClass = classOf[ContigScatterFunction]
        realigner.setupGatherFunction = {
          case (gather: BamGatherFunction, _) =>
            gather.memoryLimit = Some(6)
            gather.jarFile = qscript.picardFixMatesJar
            // Don't pass this AS=true to fix mates!
            gather.assumeSorted = None
        }
      } else {
        realigner.out = swapExt(bam,"bam","unfixed.cleaned.bam")

        // Explicitly run fix mates if the function won't be scattered.

        fixMates.memoryLimit = Some(6)
        fixMates.jarFile = qscript.picardFixMatesJar
        fixMates.unfixed = realigner.out
        fixMates.fixed = cleaned_bam
        fixMates.analysisName = "FixMates_"+bam.getName
        // Add the fix mates explicitly
      }

      var samtoolsindex = new SamtoolsIndexFunction
      samtoolsindex.bamFile = cleaned_bam
      samtoolsindex.analysisName = "index_"+cleaned_bam.getName

      if (!qscript.skip_cleaning) {
        if ( realigner.scatterCount > 1 ) {
            add(targetCreator,realigner,samtoolsindex)
        } else {
            add(targetCreator,realigner,fixMates,samtoolsindex)
        }
      }
    }

    val recalibratedBamFiles = recalibratedSamples
            .map(_.getBamFiles.get("recalibrated"))
            .toList

    val cleanBamFiles = qscript.pipeline.getSamples
            .filter(_.getBamFiles.contains("cleaned"))
            .map(_.getBamFiles.get("cleaned"))
            .toList

    // actually make calls
    if (!qscript.skip_cleaning) {
      //endToEnd(cleanedBase, cleanBamFiles, adprRscript, seq, expKind)
      endToEnd(cleanedBase, cleanBamFiles)
    } else {
      //endToEnd(uncleanedBase, recalibratedBamFiles, adprRscript, seq, expKind)
      endToEnd(uncleanedBase, recalibratedBamFiles)
    }
  }

  //def endToEnd(base: String, bamFiles: List[File], adprthing: File, seqinfo: String, exptype: String) = {
  def endToEnd(base: String, bamFiles: List[File]) = {

    // step through the un-indel-cleaned graph:
    // 1a. call snps and indels
    val snps = new UnifiedGenotyper with CommandLineGATKArgs
    snps.analysisName = base+"_SNP_calls"
    snps.input_file = bamFiles
    snps.group :+= "Standard"
    snps.out = new File(base+".vcf")
    snps.standard_min_confidence_threshold_for_emitting = Some(10)
    snps.min_mapping_quality_score = Some(20)
    snps.min_base_quality_score = Some(20)
    snps.downsample_to_coverage = Some(qscript.downsampling_coverage)
    //snps.annotation :+= "QualByDepthV2"

    if (qscript.trigger != null) {
      snps.trigger_min_confidence_threshold_for_calling = Some(30)
      snps.rodBind :+= RodBind("trigger", "VCF", qscript.trigger)
      // TODO: triggers need to get a name for comp-ing them if not dbSNP?
      snps.rodBind :+= RodBind( "compTrigger", "VCF", qscript.trigger )
    }

    // todo -- add generalize comp inputs
    //if ( qscript.comp1KGCEU != null ) {
    //  snps.rodBind :+= RodBind( "comp1KG_CEU", "VCF", qscript.comp1KGCEU )
    //}

    snps.scatterCount = qscript.num_snp_scatter_jobs

    // indel genotyper does one sample at a time
    var indelCallFiles = List.empty[RodBind]
    var indelGenotypers = List.empty[IndelGenotyperV2 with CommandLineGATKArgs]
    var loopNo = 0
    var priority = ""
    for ( bam <- bamFiles ) {
      var indel = new IndelGenotyperV2 with CommandLineGATKArgs
      indel.analysisName = "IndelGenotyper_"+bam.getName
      indel.input_file :+= bam
      indel.out = swapExt(bam,".bam",".indels.vcf")
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
    mergeIndels.out = new TaggedFile(qscript.pipeline.getProject.getName+".indels.vcf","vcf")
    mergeIndels.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    mergeIndels.priority = priority
    mergeIndels.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    mergeIndels.rodBind = indelCallFiles
    mergeIndels.analysisName = base+"_MergeIndels"

    // 1b. genomically annotate SNPs -- no longer slow
    val annotated = new GenomicAnnotator with CommandLineGATKArgs
    annotated.rodBind :+= RodBind("variant", "VCF", snps.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.refseqTable)
    //annotated.rodBind :+= RodBind("dbsnp", "AnnotatorInputTable", qscript.dbsnpTable)
    annotated.out = swapExt(snps.out,".vcf",".annotated.vcf")
    //annotated.select :+= "dbsnp.name,dbsnp.refUCSC,dbsnp.strand,dbsnp.observed,dbsnp.avHet"
    annotated.rodToIntervalTrackName = "variant"
    annotated.analysisName = base+"_GenomicAnnotator"

    // 2.a filter on cluster and near indels
    val masker = new VariantFiltration with CommandLineGATKArgs
    masker.variantVCF = annotated.out
    masker.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    masker.maskName = "NearIndel"
    masker.clusterWindowSize = Some(10)
    masker.clusterSize = Some(3)
    masker.out = swapExt(annotated.out,".vcf",".indel.masked.vcf")
    masker.analysisName = base+"_Cluster_and_Indel_filter"

    // 2.b hand filter with standard filter
    val handFilter = new VariantFiltration with CommandLineGATKArgs
    handFilter.variantVCF = masker.out
    handFilter.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    handFilter.filterName ++= List("StrandBias","AlleleBalance","QualByDepth","HomopolymerRun")
    handFilter.filterExpression ++= List("\"SB>=0.10\"","\"AB>=0.75\"","\"QD<5.0\"","\"HRun>=4\"")
    handFilter.out = swapExt(annotated.out,".vcf",".handfiltered.vcf")
    handFilter.analysisName = base+"_HandFilter"

    // 3.i generate gaussian clusters on the masked vcf
    // todo -- args for annotations?
    // todo -- args for resources (properties file)
    val clusters = new GenerateVariantClusters with CommandLineGATKArgs
    clusters.rodBind :+= RodBind("input", "VCF", masker.out)
    clusters.DBSNP = qscript.pipeline.getProject.getDbsnpFile
    val clusters_clusterFile = swapExt(new File(snps.out.getAbsolutePath),".vcf",".cluster")
    clusters.clusterFile = clusters_clusterFile
    clusters.memoryLimit = Some(6)
    clusters.jobQueue = "gsa"

    clusters.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
    clusters.analysisName = base+"_Cluster"

    // 3.ii apply gaussian clusters to the masked vcf
    val recalibrate = new VariantRecalibrator with CommandLineGATKArgs
    recalibrate.clusterFile = clusters.clusterFile
    recalibrate.DBSNP = qscript.pipeline.getProject.getDbsnpFile
    recalibrate.rodBind :+= RodBind("input", "VCF", masker.out)
    recalibrate.out = swapExt(masker.out,".vcf",".recalibrated.vcf")
    recalibrate.target_titv = qscript.target_titv
    recalibrate.report_dat_file = swapExt(masker.out,".vcf",".recalibrate.dat")
    recalibrate.tranches_file = swapExt(masker.out,".vcf",".recalibrate.tranches")
    recalibrate.analysisName = base+"_VariantRecalibrator"

    // 3.iii apply variant cuts to the clusters
    val cut = new ApplyVariantCuts with CommandLineGATKArgs
    cut.rodBind :+= RodBind("input", "VCF", recalibrate.out)
    cut.out = swapExt(recalibrate.out,".vcf",".tranched.vcf")
    cut.tranches_file = recalibrate.tranches_file
    // todo -- fdr inputs, etc
    cut.fdr_filter_level = Some(1)
    cut.analysisName = base+"_ApplyVariantCuts"

    // 4. Variant eval the cut and the hand-filtered vcf files
    val eval = new VariantEval with CommandLineGATKArgs
    eval.rodBind :+= RodBind("evalOptimized", "VCF", cut.out)
    eval.rodBind :+= RodBind("evalHandFiltered", "VCF", handFilter.out)
    eval.evalModule ++= List("CountFunctionalClasses", "CompOverlap", "CountVariants", "TiTvVariantEvaluator")
    eval.reportLocation = new File(base+".eval")
    eval.reportType = Option(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.R)
    eval.analysisName = base+"_VariantEval"

    add(snps)

    // 5. Run the ADPR and make pretty stuff

//    val adpr = new CommandLineFunction{
//     @Input(doc="Dependent files") var dependents: File = _
//     @Output(doc="Automated Data processing report") var out: File = _
//      var setname: String
//      var protocol: String
//      var sequencer: String
//      var scriptloc: File
//      def commandLine = "Rscript %s %s %s %s"
//        .format(scriptloc, setname, protocol, sequencer)
//    }
//
//    adpr.setname = base
//    adpr.scriptloc = adprthing
//    adpr.sequencer = seqinfo
//    adpr.protocol = exptype
//    adpr.dependents = eval.reportLocation
//    adpr.out = new File(base + "_adpr.pdf")
//    adpr.analysisName = base + "_ADPR"
    //In order for ADPR to finish successfully, a squid file for both the lane and sample level data needs to be
    // produced, reformatted and named <projectBase>_lanes.txt or <projectBase>_samps.txt, respectively. These files
    // to be in the working directory. When database access is ready, this and the protocol and sequencer parameters of
    //the r script will go away.

    for ( igv2 <- indelGenotypers ) {
      add(igv2)
    }

//    add(mergeIndels,annotated,masker,handFilter,clusters,recalibrate,cut,eval,adpr)
    add(mergeIndels,annotated,masker,handFilter,clusters,recalibrate,cut,eval)
  }
}
