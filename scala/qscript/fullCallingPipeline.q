import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class fullCallingPipeline extends QScript {
  qscript =>

  @Argument(doc="contigIntervals", shortName="contigIntervals")
  var contigIntervals: File = _

  @Argument(doc="numContigs", shortName="numContigs")
  var numContigs: Int = _

  @Argument(doc="project", shortName="project")
  var project: String = _

  @Input(doc="trigger", shortName="trigger", required=false)
  var trigger: File = _

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

trait CommandLineGATKArgs extends CommandLineGATK {
  this.intervals = qscript.intervals
  this.jarFile = qscript.gatkJar
}

// ------------ SETUP THE PIPELINE ----------- //

// todo -- the unclean and clean pipelines are the same, so the code can be condensed significantly

  def script = {
    val projectBase: String = qscript.project
    val cleanedBase: String = projectBase + ".cleaned"
    val uncleanedBase: String = projectBase + ".uncleaned"
// there are commands that use all the bam files
    var cleanBamFiles = List.empty[File]

for ( bam <- bamFiles ) {

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
 //realigner.targetIntervals = targetCreator.out
 realigner.targetIntervals = targetCreator.out.getAbsolutePath
 realigner.scatterCount = qscript.numContigs
 realigner.out = cleaned_bam
 realigner.scatterClass = classOf[ContigScatterFunction]
 realigner.setupGatherFunction = { case (f: BamGatherFunction, _) => f.jarFile = qscript.picardFixMatesJar }
 realigner.jobQueue = "long"

 // put clean bams in clean genotypers

 cleanBamFiles :+= realigner.out

  add(targetCreator,realigner)
}
    endToEnd(uncleanedBase,bamFiles)
    endToEnd(cleanedBase,cleanBamFiles)
  }

def endToEnd(base: String, bamFiles: List[File]) = {

 // step through the un-indel-cleaned graph:
 // 1a. call snps and indels
 val snps = new UnifiedGenotyper with CommandLineGATKArgs
 snps.input_file = bamFiles
 snps.group :+= "Standard"
 snps.annotation :+= "MyHamplotypeScore"
 snps.variants_out = new File(base+".vcf")
 snps.standard_min_confidence_threshold_for_emitting = Some(10)
 snps.min_mapping_quality_score = Some(20)
 snps.min_base_quality_score = Some(20)
 snps.downsampling_type = Some(DownsampleType.EXPERIMENTAL_BY_SAMPLE)
 snps.downsample_to_coverage = Some(200)
 // todo -- add input for comps, triggers, etc
 if (qscript.trigger != null) {
  snps.trigger_min_confidence_threshold_for_calling = Some(30)
  snps.rodBind :+= RodBind("trigger", "VCF", qscript.trigger)
 }
 // todo -- hack -- get this from the command line, or properties
 snps.rodBind :+= RodBind( "comp1KG_CEU", "VCF", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/CEU.low_coverage.2010_07.sites.hg18.vcf.gz") )


 // TODO: what is the 1KG_ALL track?
 //snps.rodBind :+= RodBind( "comp1KG_ALL", "VCF", qscript.trigger )


 snps.scatterCount = 100
 val indels = new UnifiedGenotyper with CommandLineGATKArgs
 indels.input_file = bamFiles
 indels.variants_out = new File(base+".indels.vcf")
 indels.genotype_model = Some(Model.INDELS)
 indels.scatterCount = 100
 // todo -- add inputs for the indel genotyper
 // 1b. genomically annotate SNPs -- slow, but scatter it
 val annotated = new GenomicAnnotator with CommandLineGATKArgs
 annotated.rodBind :+= RodBind("variant", "VCF", snps.variants_out)
 annotated.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.refseqTable)
 annotated.rodBind :+= RodBind("dbsnp", "AnnotatorInputTable", qscript.dbsnpTable)
 annotated.vcfOutput = swapExt(snps.variants_out,".vcf",".annotated.vcf")
 annotated.select :+= "dbsnp.name,dbsnp.refUCSC,dbsnp.strand,dbsnp.observed,dbsnp.avHet"
 annotated.rodToIntervalTrackName = "variant"
 annotated.scatterCount = 100
 // 2.a filter on cluster and near indels
 val masker = new VariantFiltration with CommandLineGATKArgs
 masker.rodBind :+= RodBind("variant", "VCF", annotated.vcfOutput)
 masker.rodBind :+= RodBind("mask", "VCF", indels.variants_out)
 masker.maskName = "NearIndel"
 masker.clusterWindowSize = Some(20)
 masker.clusterSize = Some(7)
 masker.out = swapExt(annotated.vcfOutput,".vcf",".indel.masked.vcf")
 // todo -- snp cluster args?
 // 2.b hand filter with standard filter
 val handFilter = new VariantFiltration with CommandLineGATKArgs
 handFilter.rodBind :+= RodBind("variant", "VCF", annotated.vcfOutput)
 handFilter.rodBind :+= RodBind("mask", "VCF", indels.variants_out)
 handFilter.filterName ++= List("StrandBias","AlleleBalance","QualByDepth","HomopolymerRun")
 handFilter.filterExpression ++= List("\"SB>=0.10\"","\"AB>=0.75\"","QD<5","\"HRun>=4\"")
 handFilter.out = swapExt(annotated.vcfOutput,".vcf",".handfiltered.vcf")
 // 3.i generate gaussian clusters on the masked vcf
 val clusters = new GenerateVariantClusters with CommandLineGATKArgs
 clusters.rodBind :+= RodBind("input", "VCF", masker.out)
 //clusters.clusterFile = swapExt(snps.variants_out,".vcf",".cluster")
 val clusters_clusterFile = swapExt(snps.variants_out,".vcf",".cluster")
 clusters.clusterFile = clusters_clusterFile.getAbsolutePath
 clusters.memoryLimit = Some(8)
 clusters.jobQueue = "hugemem"
 // todo -- args for annotations?
 // todo -- args for resources (properties file)
 clusters.use_annotation ++= List("QD", "SB", "MyHaplotypeScore", "HRun")
 clusters.path_to_resources = "/humgen/gsa-scr1/chartl/sting/R"
 // 3.ii apply gaussian clusters to the masked vcf
 val recalibrate = new VariantRecalibrator with CommandLineGATKArgs
 recalibrate.clusterFile = clusters.clusterFile
 recalibrate.rodBind :+= RodBind("input", "VCF", masker.out)
 recalibrate.out = swapExt(masker.out,".vcf",".optimized.vcf")
 // todo -- inputs for Ti/Tv expectation and other things
 recalibrate.target_titv = Some(2.1)
 // 3.iii apply variant cuts to the clusters
 val cut = new ApplyVariantCuts with CommandLineGATKArgs
 cut.rodBind :+= RodBind("input", "VCF", recalibrate.out)
 //cut.outputVCFFile = swapExt(recalibrate.out,".vcf",".tranched.vcf")
 //cut.tranchesFile = swapExt(recalibrate.out,".vcf",".tranch")
 val cut_outputVCFFile = swapExt(recalibrate.out,".vcf",".tranched.vcf")
 val cut_tranchesFile = swapExt(recalibrate.out,".vcf",".tranch")
 cut.outputVCFFile = cut_outputVCFFile.getAbsolutePath
 cut.tranchesFile = cut_tranchesFile.getAbsolutePath
 // todo -- fdr inputs, etc
 cut.fdr_filter_level = Some(10)
 // 4. Variant eval the cut and the hand-filtered vcf files
 val eval = new VariantEval with CommandLineGATKArgs
 eval.rodBind :+= RodBind("evalOptimized", "VCF", cut_outputVCFFile)
 eval.rodBind :+= RodBind("evalHandFiltered", "VCF", handFilter.out)
 // todo -- make comp tracks command-line arguments or properties
 eval.evalModule ++= List("CountFunctionalClasses", "CompOverlap", "CountVariants", "TiTvVariantEvaluator")
 eval.out = new File(base+".eval")

 add(snps,indels,annotated,masker,handFilter,clusters,recalibrate,cut,eval)
}

}
