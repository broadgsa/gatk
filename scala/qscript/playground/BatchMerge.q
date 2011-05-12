import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.walkers.genotyper.{GenotypeLikelihoodsCalculationModel, UnifiedGenotyperEngine}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.vcf.{VCFSimpleMerge, VCFExtractSites,VCFExtractIntervals}
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.utils.text.XReadLines

class batchMergePipeline extends QScript {
  batchMerge =>

  @Argument(doc="VCF list",shortName="vcfs") var vcfList: File = _
  @Argument(doc="bam list",shortName="bams") var bamList: File = _
  @Argument(doc="sting dir",shortName="sting") var stingDir: String = _
  @Argument(doc="reference file",shortName="ref") var ref: File = _
  @Argument(doc="batched output",shortName="batch") var batchOut: File = _
  //@Argument(doc="read UG settings from header",shortName="ugh") var ugSettingsFromHeader : Boolean = false
  @Hidden @Argument(doc="Min base q",shortName="mbq",required=false) var mbq : Int = 20
  @Hidden @Argument(doc="Min map q",shortName="mmq",required=false) var mmq : Int = 20
  @Hidden @Argument(doc="baq gap open penalty, using sets baq to calc when necessary",shortName="baqp",required=false) var baq : Int = -1

  def script = {

    var vcfs : List[File] = extractFileEntries(vcfList)
    var bams : List[File] = extractFileEntries(bamList)

    trait ExtractArgs extends VCFExtractSites {
      this.keepFilters = false
      this.keepInfo = false
      this.keepQual = false
    }



    trait CombineVariantsArgs extends CombineVariants {
      this.reference_sequence = batchMerge.ref
      this.jarFile = new File(batchMerge.stingDir+"/dist/GenomeAnalysisTK.jar")
      this.scatterCount = 10
      this.memoryLimit=4
    }
    var combine : CombineVariants = new CombineVariants with CombineVariantsArgs
    combine.out = swapExt(batchOut,".vcf",".variant.combined.vcf")
    combine.rodBind ++= vcfs.map( u => new RodBind(u.getName,"vcf",u) )

    add(combine)

    var getVariantAlleles : List[VCFExtractSites] = vcfs.map( u => new VCFExtractSites(u, swapExt(batchOut.getParent,u,".vcf",".alleles.vcf")) with ExtractArgs)
    var combineVCFs : VCFSimpleMerge = new VCFSimpleMerge
    combineVCFs.vcfs = getVariantAlleles.map(u => u.outVCF)
    combineVCFs.fai = new File(ref.getAbsolutePath+".fai")
    combineVCFs.outVCF = swapExt(batchOut,".vcf",".pf.alleles.vcf")
    var extractIntervals : VCFExtractIntervals = new VCFExtractIntervals(combine.out,swapExt(combine.out,".vcf",".intervals.list"),true)
    //addAll(getVariantAlleles)
    //add(combineVCFs,extractIntervals)
    add(extractIntervals)

    trait CalcLikelihoodArgs extends UGCalcLikelihoods {
      this.reference_sequence = batchMerge.ref
      this.min_base_quality_score = batchMerge.mbq
      this.min_mapping_quality_score = batchMerge.mmq
      if ( batchMerge.baq >= 0 ) {
        this.baqGapOpenPenalty = batchMerge.baq
        this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
      }
      this.intervals :+= extractIntervals.listOut
      this.allelesVCF = combine.out
      this.jarFile = new File(stingDir+"/dist/GenomeAnalysisTK.jar")
      this.memoryLimit = 4
      this.scatterCount = 60
      this.output_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
      this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    }

    def newUGCL( bams: (List[File],Int) ) : UGCalcLikelihoods = {
      var ugcl = new UGCalcLikelihoods with CalcLikelihoodArgs
      ugcl.input_file ++= bams._1
      ugcl.out = new File("MBatch%d.likelihoods.vcf".format(bams._2))
      return ugcl
    }

    var calcs: List[UGCalcLikelihoods] = bams.grouped(20).toList.zipWithIndex.map(u => newUGCL(u))
    addAll(calcs)

    trait CallVariantsArgs extends UGCallVariants {
      this.reference_sequence = batchMerge.ref
      this.intervals :+= extractIntervals.listOut
      this.jarFile = new File(stingDir+"/dist/GenomeAnalysisTK.jar")
      this.scatterCount = 30
      this.memoryLimit = 8
      this.output_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
      this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    }

    var cVars : UGCallVariants = new UGCallVariants with CallVariantsArgs
    cVars.rodBind ++= calcs.map( a => new RodBind("variant"+a.out.getName.replace(".vcf",""),"vcf",a.out) )
    cVars.out = batchOut
    add(cVars)
  }

  override def extractFileEntries(in: File): List[File] = {
    return (new XReadLines(in)).readLines.toList.map( new File(_) )
  }
}
