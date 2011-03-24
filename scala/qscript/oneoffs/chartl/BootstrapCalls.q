import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.{GenotypeMergeType, VariantMergeType}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class BootstrapCalls extends QScript {
  @Argument(doc="Bam file list",shortName="I",required=true)
  var bamList: File = _
  @Argument(doc="Intervals file",shortName="L",required=true)
  var intervalFile: File = _
  @Argument(doc="Output file",shortName="o",required=true)
  var bootstrapMergedOut: File = _
  @Argument(doc="Reference file",shortName="R",required=true)
  var reference: File = _
  @Argument(doc="Downsampling Level",shortName="D",required=false)
  var downsamplingLevel: Int = 4
  @Argument(doc="Num Bootstrap Callsets",shortName="B",required=false)
  var numberOfBootstraps: Int = 25
  @Argument(doc="call confidence",shortName="conf",required=false)
  var standCallConf: Double = 4.0
  @Argument(doc="dbsnp file (vcf version)",shortName="dbsnp",required=false)
  var dbsnp: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  @Argument(doc="sting jar",shortName="s",required=true)
  var sting: File = _

  /**********************
   *            URGENT NOTE:
   *            for this to do any good you need to take out the random seeds in
   *            ReservoirDownsampler: 20
   *            MathUtils: 649
   *
   *            You will also need to hack the recalibrator to always trust AC (which are no longer integer-valued)
   *            and to deal with double-valued AC fields
   */

  def script = {
    val bams: List[File] = extractFileEntries(bamList)
    trait UGArgs extends UnifiedGenotyper {
      this.input_file = bams
      this.reference_sequence = reference
      this.dcov = downsamplingLevel
      this.intervals :+= intervalFile
      this.stand_call_conf = standCallConf
      this.stand_emit_conf = standCallConf
      this.rodBind :+= new RodBind("dbsnp","vcf",dbsnp)
      this.scatterCount = 20
      this.jarFile = sting
      this.memoryLimit = 4
    }

    val bootstrapBase = swapExt(bootstrapMergedOut,".vcf",".boot%d.vcf").getAbsolutePath
    var calls : List[UnifiedGenotyper] = Nil
    for ( i <- 0 until (numberOfBootstraps+1) ) {
      var ug : UnifiedGenotyper = new UnifiedGenotyper with UGArgs
      ug.out = new File(bootstrapBase.format(i))
      ug.analysisName = "Boostrap%d".format(i)
      calls :+= ug
    }

    addAll(calls)

    trait MergeArgs extends BootstrapCallsMerger {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.scatterCount = 40
      this.jarFile = sting
      this.memoryLimit = 4
      this.rodBind ++= calls.map(u => u.out).zipWithIndex.map(u => new RodBind("bootstrap_%d".format(u._2),"vcf",u._1))
      this.out = bootstrapMergedOut
    }

    var merge : BootstrapCallsMerger = new BootstrapCallsMerger with MergeArgs
    add(merge)

    trait ClusterArgs extends GenerateVariantClusters {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.rodBind :+= new RodBind("input","vcf",merge.out)
      this.rodBind :+= new RodBind("hapmap","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"))
      this.rodBind :+= new RodBind("truthHapMap","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"))
      this.rodBind :+= new RodBind("1kg","vcf", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.sites.vcf"))
      this.rodBind :+= new RodBind("truth1kg","vcf", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.sites.vcf"))
      this.cluster_file = swapExt(bootstrapMergedOut,"vcf","cluster")
      this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      this.qual = 100
      this.std = 3.5
      this.mG = 8
      this.trustAllPolymorphic = true
      this.memoryLimit = 8
      this.jarFile = sting
    }

    var clust : GenerateVariantClusters = new GenerateVariantClusters with ClusterArgs
    add(clust)

    trait VQSRArgs extends VariantRecalibrator {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.out = swapExt(bootstrapMergedOut,"vcf","recal.vcf")
      this.rodBind :+= new RodBind("input","vcf",merge.out)
      this.rodBind :+= new RodBind("hapmap","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"))
      this.rodBind :+= new RodBind("1kg","vcf", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.sites.vcf"))
      this.rodBind :+= new RodBind("truthHapMap","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"))
      this.rodBind :+= new RodBind("truth1kg","vcf", new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/1212samples.b37.sites.vcf"))
      this.cluster_file = swapExt(bootstrapMergedOut,"vcf","cluster")
      this.sm = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY
      this.tranche ++= List("0.1", "0.5", "0.7", "1.0", "3.0", "5.0", "10.0", "100.0")
      this.trustAllPolymorphic = true
      this.tranchesFile = swapExt(bootstrapMergedOut,"vcf","tranche")
      this.memoryLimit=8
      this.jarFile = sting
      this.rodBind :+= new RodBind("dbsnp","vcf",dbsnp)
    }

    var recal : VariantRecalibrator = new VariantRecalibrator with VQSRArgs
    add(recal)

    trait CutArgs extends ApplyVariantCuts {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.rodBind :+= new RodBind("input","vcf",recal.out)
      this.tranchesFile = recal.tranchesFile
      this.fdr_filter_level = 1.0
      this.out = swapExt(bootstrapMergedOut,".vcf",".recal.cut.vcf")
      this.jarFile = sting
      this.memoryLimit = 4
      this.scatterCount = 5
    }

    var cut : ApplyVariantCuts = new ApplyVariantCuts with CutArgs
    add(cut)

    class RmHeader extends CommandLineFunction {
      @Input(doc="vcf") var vcf : File = _
      @Output(doc="headerless vcf") var noheadvcf : File = _

      def commandLine : String = {
        "head -n 1 %s > %s ; grep -v \\#\\# %s >> %s".format(vcf.getAbsolutePath,noheadvcf.getAbsolutePath,
          vcf.getAbsolutePath,noheadvcf.getAbsolutePath)
      }
    }

    var rm : RmHeader = new RmHeader
    rm.vcf = cut.out
    rm.noheadvcf = swapExt(cut.out,".vcf",".nohead.vcf")
    add(rm)

    trait CombineArgs extends CombineVariants {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.rodBind :+= new RodBind("loCov","vcf",rm.noheadvcf)
      this.rodBind :+= new RodBind("hiCov","vcf",new File("/humgen/gsa-pipeline/PVQF4/all_batches_v001/batch_001/SnpCalls/ESPGO_Gabriel_NHLBI_EOMI_setone_EOMI_Project.cleaned.annotated.handfiltered.vcf"))
      this.variantMergeOptions = VariantMergeType.UNION
      this.genotypeMergeOptions = GenotypeMergeType.PRIORITIZE
      this.priority = "hiCov,loCov"
      this.out = swapExt(bootstrapMergedOut,".vcf",".merged.combined.vcf")
      this.jarFile = sting
      this.memoryLimit = 6
    }

    var combine : CombineVariants = new CombineVariants with CombineArgs
    add(combine)

    trait EvalArgs extends VariantEval {
      this.reference_sequence = reference
      this.intervals :+= intervalFile
      this.rodBind :+= new RodBind("evalCombined","vcf",combine.out)
      //this.rodBind :+= new RodBind("evalCut","vcf",rm.noheadvcf)
      //this.rodBind :+= new RodBind("evalFCP","vcf",new File("/humgen/gsa-pipeline/PVQF4/all_batches_v001/batch_001/SnpCalls/ESPGO_Gabriel_NHLBI_EOMI_setone_EOMI_Project.cleaned.annotated.handfiltered.vcf"))
      this.rodBind :+= new RodBind("dbsnp","vcf",dbsnp)
      this.jarFile = sting
      this.ST = List("Filter","Novelty","JexlExpression")
      this.select_names = List("lowOnly","filteredInLow","Intersection","filteredInHi","hiOnly","filteredInAll")
      this.select_exps = List("\"set == 'loCov'\"","\"set == 'hiCov-filterInloCov'\"",
                              "\"set == 'Intersection'\"", "\"set == 'filterInhiCov-loCov'\"",
                              "\"set == 'hiCov'\"","\"set == 'FilteredInAll'\"")
      this.EV = List("TiTvVariantEvaluator","CountVariants","CompOverlap")
      this.out = swapExt(bootstrapMergedOut,".vcf",".merged.combined.eval")
      this.nt = 8
      this.memoryLimit = 12
    }

    var eval : VariantEval = new VariantEval with EvalArgs
    add(eval)
  }
}
