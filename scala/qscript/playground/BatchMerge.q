import java.io.{FileReader, BufferedReader}
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.vcf.{VCFSimpleMerge, VCFExtractSites,VCFExtractIntervals}
import org.broadinstitute.sting.queue.pipeline.{ProjectManagement, BamProcessing, VariantCalling}
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.utils.yaml.YamlUtils

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
  @Hidden @Argument(doc="Max mismatching bases",shortName="mmb",required=false) var mmb : Int = 3
  @Hidden @Argument(doc="baq gap open penalty, using sets baq to calc when necessary",shortName="baqp",required=false) var baq : Int = -1

  def script = {

    var vcfs : List[File] = extractFileEntries(vcfList)
    var bams : List[File] = extractFileEntries(bamList)

    trait ExtractArgs extends VCFExtractSites {
      this.keepFilters = true
      this.keepInfo = false
      this.keepQual = false
    }

    var getVariantAlleles : List[VCFExtractSites] = vcfs.map( u => new VCFExtractSites(u, swapExt(batchOut.getParent,u,".vcf",".alleles.vcf")) with ExtractArgs)
    var combineVCFs : VCFSimpleMerge = new VCFSimpleMerge
    combineVCFs.vcfs = getVariantAlleles.map(u => u.outVCF)
    combineVCFs.fai = new File(ref.getAbsolutePath+".fai")
    combineVCFs.outVCF = swapExt(batchOut,".vcf",".pf.alleles.vcf")
    var extractIntervals : VCFExtractIntervals = new VCFExtractIntervals(combineVCFs.outVCF,swapExt(combineVCFs.outVCF,".vcf",".intervals.list"),true)
    addAll(getVariantAlleles)
    add(combineVCFs,extractIntervals)

    trait CalcLikelihoodArgs extends UGCalcLikelihoods {
      this.reference_sequence = batchMerge.ref
      this.max_mismatches_in_40bp_window = Some(batchMerge.mmb)
      this.min_base_quality_score = Some(batchMerge.mbq)
      this.min_mapping_quality_score = Some(batchMerge.mmq)
      if ( batchMerge.baq >= 0 ) {
        this.baqGapOpenPenalty = Some(batchMerge.baq)
        this.baq = Some(BAQ.CalculationMode.CALCULATE_AS_NECESSARY)
      }
      this.intervals :+= extractIntervals.listOut
      this.alleleVCF = combineVCFs.outVCF
      this.output_all_callable_bases = true
      this.jarFile = new File(stingDir+"/dist/GenomeAnalysisTK.jar")
      this.memoryLimit = Some(4)
      this.scatterCount = 60
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
      this.output_all_callable_bases = true
      this.reference_sequence = batchMerge.ref
      this.intervals :+= extractIntervals.listOut
      this.jarFile = new File(stingDir+"/dist/GenomeAnalysisTK.jar")
      this.scatterCount = 30
      this.memoryLimit = Some(8)
    }

    var cVars : UGCallVariants = new UGCallVariants with CallVariantsArgs
    cVars.rodBind ++= calcs.map( a => new RodBind("variant"+a.out.getName.replace(".vcf",""),"vcf",a.out) )
    cVars.out = batchOut
    add(cVars)
  }

  def extractFileEntries(in: File): List[File] = {
    return (new XReadLines(in)).readLines.toList.map( new File(_) )
  }
}
