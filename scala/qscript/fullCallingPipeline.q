import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.QScript

class fullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc = "reference", shortName="R")
  var reference: File = _

  @Argument(doc="contigIntervals", shortName="contigIntervals")
  var contigIntervals: File = _

  @Argument(doc="numContigs", shortName="numContigs")
  var numContigs: Int = _

  @Argument(doc="project", shortName="project")
  var project: String = _

  @Input(doc="trigger", shortName="trigger", required=false)
  var trigger: File = _

  @Input(doc="compCEU",shortName="ceu",required=false)
  var comp1KGCEU: File = _

  @Input(doc="refseqTable", shortName="refseqTable")
  var refseqTable: File = _

  @Input(doc="dbsnpTable", shortName="dbsnpTable")
  var dbsnpTable: File = _

  @Input(doc="Picard FixMateInformation.jar.  At the Broad this can be found at /seq/software/picard/current/bin/FixMateInformation.jar.  Outside the broad see http://picard.sourceforge.net/")
  var picardFixMatesJar: File = _

  @Input(doc="intervals")
  var intervals: File = _

  @Input(doc="bam files", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="gatk jar")
  var gatkJar: File = _

  @Input(doc="SNP cluster filter -- number SNPs",shortName="snpsInCluster",required=false)
  var snpsInCluster = 4

  @Input(doc="SNP cluster filter -- window size",shortName="snpClusterWindow",required=false)
  var snpClusterWindow = 7

  @Input(doc="dbSNP version",shortName="D")
  var dbSNP: File = _

  @Input(doc="target titv for recalibration",shortName="titv",required=false)
  var target_titv = 2.1

  @Input(doc="downsampling coverage",shortName="dcov",required=false)
  var downsampling_coverage = 200

  @Input(doc="Number of jobs to scatter unifeid genotyper",shortName="snpScatter",required=false)
  var num_snp_scatter_jobs = 50

  @Input(doc="Number of jobs to scatter indel genotyper",shortName="indelScatter",required=false)
  var num_indel_scatter_jobs = 5


  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals = qscript.intervals
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
  }


  // ------------ SETUP THE PIPELINE ----------- //


  def script = {
    val projectBase: String = qscript.project
    val cleanedBase: String = projectBase + ".cleaned"
    val uncleanedBase: String = projectBase + ".uncleaned"
    // there are commands that use all the bam files
    var cleanBamFiles = List.empty[File]

    for ( bam <- qscript.bamFiles ) {

      // put unclean bams in unclean genotypers

      // in advance, create the extension files

      val indel_targets = swapExt(bam,"bam","realigner_targets.interval_list")
      val cleaned_bam = swapExt(bam,"bam","cleaned.bam") // note-- the scatter is in the definition itself

      // create the cleaning commands

      val targetCreator = new RealignerTargetCreator with CommandLineGATKArgs
      targetCreator.input_file :+= bam
      targetCreator.out = indel_targets

      val realigner = new IndelRealigner with CommandLineGATKArgs
      realigner.input_file = targetCreator.input_file
      realigner.intervals = qscript.contigIntervals
      realigner.targetIntervals = new java.io.File(targetCreator.out.getAbsolutePath)
      realigner.scatterCount = qscript.numContigs
      realigner.out = cleaned_bam
      realigner.scatterClass = classOf[ContigScatterFunction]
      realigner.setupGatherFunction = { case (f: BamGatherFunction, _) => f.jarFile = qscript.picardFixMatesJar }
      realigner.jobQueue = "week"

      var samtoolsindex = new SamtoolsIndexFunction
      samtoolsindex.bamFile = cleaned_bam

      // put clean bams in clean genotypers

      cleanBamFiles :+= realigner.out

      add(targetCreator,realigner,samtoolsindex)
    }

    // actually make calls
    endToEnd(uncleanedBase,qscript.bamFiles)
    endToEnd(cleanedBase,cleanBamFiles)
  }

  def endToEnd(base: String, bamFiles: List[File]) = {

    // step through the un-indel-cleaned graph:
    // 1a. call snps and indels
    val snps = new UnifiedGenotyper with CommandLineGATKArgs
    snps.input_file = bamFiles
    snps.group :+= "Standard"
    snps.out = new File(base+".vcf")
    snps.standard_min_confidence_threshold_for_emitting = Some(10)
    snps.min_mapping_quality_score = Some(20)
    snps.min_base_quality_score = Some(20)
    snps.downsample_to_coverage = Some(200)
    snps.annotation :+= "QualByDepthV2"

    if (qscript.trigger != null) {
      snps.trigger_min_confidence_threshold_for_calling = Some(30)
      snps.rodBind :+= RodBind("trigger", "VCF", qscript.trigger)
      // TODO: triggers need to get a name for comp-ing them if not dbSNP?
      snps.rodBind :+= RodBind( "compTrigger", "VCF", qscript.trigger )
    }

    // todo -- add generalize comp inputs
    if ( qscript.comp1KGCEU != null ) {
      snps.rodBind :+= RodBind( "comp1KG_CEU", "VCF", qscript.comp1KGCEU )
    }

    snps.scatterCount = qscript.num_snp_scatter_jobs


    // indel genotyper does one sample at a time
    var indelCallFiles = List.empty[RodBind]
    var indelGenotypers = List.empty[IndelGenotyperV2 with CommandLineGATKArgs]
    var loopNo = 0
    var priority = ""
    for ( bam <- bamFiles ) {
      var indel = new IndelGenotyperV2 with CommandLineGATKArgs
      indel.input_file :+= bam
      indel.out = swapExt(bam,".bam",".indels.vcf")
      indel.downsample_to_coverage = Some(500)
      indelCallFiles :+= RodBind("v"+loopNo.toString, "VCF", indel.out)
      indel.scatterCount = qscript.num_indel_scatter_jobs

      indelGenotypers :+= indel

      if ( loopNo == 0 ) {
        priority = "v0"
      } else {
        priority += ",v"+loopNo.toString
      }
      loopNo += 1
    }
    val mergeIndels = new CombineVariants with CommandLineGATKArgs
    mergeIndels.out = new File(qscript.project+".indels.vcf")
    mergeIndels.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE)
    mergeIndels.priority = priority
    mergeIndels.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    mergeIndels.rodBind = indelCallFiles


    // 1b. genomically annotate SNPs -- no longer slow
    val annotated = new GenomicAnnotator with CommandLineGATKArgs
    annotated.rodBind :+= RodBind("variant", "VCF", snps.out)
    annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.refseqTable)
    annotated.rodBind :+= RodBind("dbsnp", "AnnotatorInputTable", qscript.dbsnpTable)
    annotated.out = swapExt(snps.out,".vcf",".annotated.vcf")
    annotated.select :+= "dbsnp.name,dbsnp.refUCSC,dbsnp.strand,dbsnp.observed,dbsnp.avHet"
    annotated.rodToIntervalTrackName = "variant"


    // 2.a filter on cluster and near indels
    val masker = new VariantFiltration with CommandLineGATKArgs
    masker.rodBind :+= RodBind("variant", "VCF", annotated.out)
    masker.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    masker.maskName = "NearIndel"
    masker.clusterWindowSize = Some(qscript.snpClusterWindow)
    masker.clusterSize = Some(qscript.snpsInCluster)
    masker.out = swapExt(annotated.out,".vcf",".indel.masked.vcf")


    // 2.b hand filter with standard filter
    val handFilter = new VariantFiltration with CommandLineGATKArgs
    handFilter.rodBind :+= RodBind("variant", "VCF", annotated.out)
    handFilter.rodBind :+= RodBind("mask", "VCF", mergeIndels.out)
    handFilter.filterName ++= List("StrandBias","AlleleBalance","QualByDepth","HomopolymerRun")
    handFilter.filterExpression ++= List("\"SB>=0.10\"","\"AB>=0.75\"","QD<5","\"HRun>=4\"")
    handFilter.out = swapExt(annotated.out,".vcf",".handfiltered.vcf")


    // 3.i generate gaussian clusters on the masked vcf
    // todo -- args for annotations?
    // todo -- args for resources (properties file)
    val clusters = new GenerateVariantClusters with CommandLineGATKArgs
    clusters.rodBind :+= RodBind("input", "VCF", masker.out)
    val clusters_clusterFile = swapExt(new File(snps.out.getAbsolutePath),".vcf",".cluster")
    clusters.clusterFile = clusters_clusterFile
    clusters.memoryLimit = Some(8)
    clusters.jobQueue = "hugemem"

    clusters.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")


    // 3.ii apply gaussian clusters to the masked vcf
    val recalibrate = new VariantRecalibrator with CommandLineGATKArgs
    recalibrate.clusterFile = clusters.clusterFile
    recalibrate.rodBind :+= RodBind("input", "VCF", masker.out)
    recalibrate.out = swapExt(masker.out,".vcf",".recalibrated.vcf")
    recalibrate.target_titv = qscript.target_titv
    recalibrate.report_dat_file = swapExt(masker.out,".vcf",".recalibrate.dat")
    recalibrate.tranches_file = swapExt(masker.out,".vcf",".recalibrate.tranches")

    // 3.iii apply variant cuts to the clusters
    val cut = new ApplyVariantCuts with CommandLineGATKArgs
    cut.rodBind :+= RodBind("input", "VCF", recalibrate.out)
    cut.out = swapExt(recalibrate.out,".vcf",".tranched.vcf")
    cut.tranches_file = recalibrate.tranches_file
    // todo -- fdr inputs, etc
    cut.fdr_filter_level = Some(1)

    
    // 4. Variant eval the cut and the hand-filtered vcf files
    val eval = new VariantEval with CommandLineGATKArgs
    eval.rodBind :+= RodBind("evalOptimized", "VCF", cut.out)
    eval.rodBind :+= RodBind("evalHandFiltered", "VCF", handFilter.out)
    eval.evalModule ++= List("CountFunctionalClasses", "CompOverlap", "CountVariants", "TiTvVariantEvaluator")
    eval.out = new File(base+".eval")

    add(snps)

    for ( igv2 <- indelGenotypers ) {
      add(igv2)
    }

    add(mergeIndels,annotated,masker,handFilter,clusters,recalibrate,cut,eval)

  }

}
