package org.broadinstitute.sting.queue.pipeline

import java.io.File
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
class VariantCalling(yaml: File,gatkJar: File) {
  vc =>

  // load attributes
  var attributes: Pipeline = YamlUtils.load(classOf[Pipeline], yaml)
  
  /**
   * Trait to propagate basic attributes throughout the library
   */
  trait StandardCommandLineGATK extends CommandLineGATK {
    this.reference_sequence = vc.attributes.getProject.getReferenceFile
    this.intervals = vc.attributes.getProject.getIntervalList
    this.DBSNP = vc.attributes.getProject.getDbsnpFile
    // set global memory limit on the low side. Additional input bams will affect it.
    this.memoryLimit = Some(2)
    this.jarFile = vc.gatkJar
  }

  /**
   * @Doc: Creates a standard UnifiedGenotyper CLF from input bams and an output file
   * @Return: UnifiedGenotyper with the standard GSA arguments
   * @TODO: Add a formula: f(#bams)=memory; allow yaml to specify triggers and perhaps other information
   */
  def StandardUnifiedGenotyper(bams : List[File], output : File) : UnifiedGenotyper = {
    var ug = new UnifiedGenotyper with StandardCommandLineGATK
    ug.analysisName = "UnifiedGenotyper"
    ug.input_file = bams
    ug.group :+= "Standard"
    ug.out = output
    ug.min_base_quality_score = Some(10)
    ug.min_mapping_quality_score = Some(10)
    ug.cap_base_quality_by_mapping_quality = true
    ug.standard_min_confidence_threshold_for_emitting = Some(10)
    ug.standard_min_confidence_threshold_for_calling = Some(30)
    ug.trigger_min_confidence_threshold_for_calling = Some(0)
    ug.downsample_to_coverage = Some(300)
    ug.dt = Some(DownsampleType.BY_SAMPLE)
    ug.scatterCount = 50

    if ( bams.size > 125 ) {
      ug.memoryLimit = Some(4)
    }

    return ug

  }

  /**
   * @Doc: Creates a CLF to call indels on a specific .bam file, outputting to a given output file
   * @Returns: An IndelGenotyperV2 CLF with standard GSA arguments
   */
  def StandardIndelGenotyper(bam : File, output: File) : IndelGenotyperV2 = {
    var ig = new IndelGenotyperV2 with StandardCommandLineGATK
    ig.analysisName = "IndelGenotyper"
    ig.input_file :+= bam
    ig.out = output
    ig.downsample_to_coverage = Some(300)
    return ig
  }

  /**
   * @Doc: Accessor method to StandardIndelGenotyper that allows it to be marked as a scatter job in a pipeline
   * @Returns: An IndelGenotyperV2 CLF with standard GSA arguments, marked as a scatter job
   */
  private def StandardIndelGenotyper(bam : File, output: File, gather: Boolean) : IndelGenotyperV2 = {
    var ig = StandardIndelGenotyper(bam,output)
    ig.isGather = gather
    return ig
  }

  /**
   * @Doc: Combines a list of indel VCFs to a single output file
   * @Returns: A CombineVariants CLF with standard GSA arguments
   */
  def StandardIndelCombine( igList : List[IndelGenotyperV2], output : File ) : CombineVariants = {
    var cv = new CombineVariants with StandardCommandLineGATK
    cv.out = output
    cv.genotypemergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    cv.variantmergeoption = Some(org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION)
    cv.analysisName = "IndelGenotyper"
    cv.isGather = true
    cv.priority = (igList.map[String,List[String]](ig => swapExt(ig.out,".vcf","").getAbsolutePath)).mkString(",")
    //cv.priority = (igList.foldLeft[List[String]](Nil)( (prLs, ig) =>  prLs ::: List(swapExt(ig.out,".vcf","").getAbsolutePath))).mkString(",")
    cv.rodBind = igList.map[RodBind,List[RodBind]](ig => new RodBind(swapExt(ig.out,".vcf","").getName,"VCF",ig.out))

    return cv
    
  }

  /**
   * @Doc: Generates indel calls on a list of .bam files, and merges those calls into an output file. This is a small pipeline.
   * @Returns: A list of CLFs that run indel calls and indel merging. User has zero control over individual indel VCF names.
   */
  def StandardIndelCalls ( bams : List[File], output : File ) : List[CommandLineGATK] = {
    var genotypers = bams.foldLeft[List[IndelGenotyperV2]](Nil)( (igs,bam) => igs ::: List(this.StandardIndelGenotyper(bam, swapExt(bam,".bam",".indels.vcf"), false)))
    var combine = this.StandardIndelCombine( genotypers, output )
    var callFunctions: List[CommandLineGATK] = genotypers
    callFunctions = combine :: callFunctions

    return callFunctions    
  }

  def StandardFilterAtIndels ( snps: File, indels: File, output : File ) : VariantFiltration = {
    var iFil = new VariantFiltration with StandardCommandLineGATK
    iFil.analysisName = "FilterAtIndels"
    iFil.out = output
    iFil.mask = indels.getAbsolutePath
    iFil.maskName = "Indel_Mask"
    iFil.variantVCF = snps
    // todo -- cluster size varies with # bams
    iFil.clusterSize = Some(5)
    iFil.clusterWindowSize = Some(8)

    return iFil
  }

  def StandardHandfilter( snps: File, output: File ) : VariantFiltration = {
    var hFil = new VariantFiltration with StandardCommandLineGATK
    hFil.analysisName = "HandFilter"
    hFil.out = output
    hFil.variantVCF = snps
    hFil.filterExpression :+= "QD<5"
    hFil.filterName :+= "LowQualByDepth"
    hFil.filterExpression :+= "AB>0.75"
    hFil.filterName :+= "HighAlleleBalance"
    hFil.filterExpression :+= "SB>-0.10"
    hFil.filterName :+= "HighStrandBias"
    
    return hFil
  }

  def StandardVariantCluster( snps: File, output: File ) : GenerateVariantClusters = {
    var genC = new GenerateVariantClusters with StandardCommandLineGATK
    genC.analysisName = "VariantQualityRecalibration"
    genC.rodBind :+= new RodBind("input","VCF",snps)
    genC.cluster_file = output
    genC.use_annotation :+= "QD"
    genC.use_annotation :+= "SB"
    genC.use_annotation :+= "HaplotypeScore"
    genC.use_annotation :+= "AB"
    genC.use_annotation :+= "HRun"

    return genC
  }

  def StandardVariantRecalibrator ( raw_vcf: File, cluster: File, target_titv: scala.Double, out_vcf: File,
                                              out_tranches: File, out_dat: File) : VariantRecalibrator = {
    var vr = new VariantRecalibrator with StandardCommandLineGATK
    vr.analysisName = "VariantQualityRecalibration"
    vr.rodBind :+= new RodBind("input","VCF",raw_vcf)
    vr.cluster_file = cluster
    vr.target_titv = target_titv
    vr.out = out_vcf
    vr.tranches_file = out_tranches
    vr.report_dat_file = out_dat
    vr.tranche :+= "0.1"
    vr.tranche :+= "1"
    vr.tranche :+= "5"
    vr.tranche :+= "10"

    return vr

  }

  def StandardApplyVariantCuts( snpRecal: File, tranches: File, output: File) : ApplyVariantCuts = {
    var avc = new ApplyVariantCuts with StandardCommandLineGATK
    avc.analysisName = "VariantQualityRecalibration"
    avc.rodBind :+= new RodBind("input","VCF",snpRecal)
    avc.out = output
    avc.tranches_file = tranches
    avc.fdr_filter_level = Some(5)
  
    return avc
  }

  def StandardRecalibrateVariants( snps: File, targetTiTv: scala.Double, recalVcf: File) : List[CommandLineGATK] = {
    var clust = StandardVariantCluster(snps, swapExt(snps,".vcf",".cluster"))
    var recal = StandardVariantRecalibrator(snps,clust.clusterFile,targetTiTv,swapExt(snps,".vcf",".recal.vcf"),
                swapExt(snps,".vcf",".recal.tranch"),swapExt(snps,".vcf",".recal.dat"))
    var cut = StandardApplyVariantCuts(recal.out,recal.tranches_file,recal.out)

    var cmds: List[CommandLineGATK] = Nil
    cmds :+= clust
    cmds :+= recal
    cmds :+= cut

    return cmds
  }

  def StandardGenomicAnnotation ( snps: File, refseqFile: File, outputVCF: File) : GenomicAnnotator = {
    var ga = new GenomicAnnotator with StandardCommandLineGATK
    ga.analysisName = "GenomicAnnotator"
    ga.variantVCF = snps
    ga.rodBind :+= new RodBind("refseq","AnnotatorInputTable",refseqFile)
    ga.rodToIntervalTrackName = "variant"
    ga.out = outputVCF

    return ga
  }

  def StandardSNPCalls( bams: List[File], output: File, targetTiTv: scala.Double, refGene: File = null ) : List[CommandLineGATK] = {
    var commands : List[CommandLineGATK] = Nil
    var dir = ""

    if ( output.getParent != null ) {
      dir = output.getParent+"/"
    }

    var raw_snp = new File(dir+vc.attributes.getProject+".raw_snps.vcf")
    var ug = StandardUnifiedGenotyper(bams, raw_snp)

    commands :+= ug

    var raw_indel = new File(dir+vc.attributes.getProject+".raw_indels.vcf")
    var ig = StandardIndelCalls(bams,raw_indel)

    commands ++= ig

    var prefilt_snp = swapExt(raw_snp,".vcf",".indel_filtered.vcf")
    var iFilt = StandardFilterAtIndels(raw_snp,raw_indel,prefilt_snp)

    commands :+= iFilt

    var annoSNP = prefilt_snp
    if ( refGene != null ) {
      annoSNP = swapExt(prefilt_snp,".vcf",".annotated.vcf")
      var annotate = StandardGenomicAnnotation(prefilt_snp, refGene, annoSNP)

      commands :+= annotate
    }

    var recal = StandardRecalibrateVariants(annoSNP, targetTiTv, output)

    commands ++= recal

    return commands
  }

  def StandardSNPCallsBothFilterTypes(bams: List[File], recalOut: File, handFilteredOut: File, targetTiTv: scala.Double, refGene: File = null ) : List[CommandLineGATK] = {
    var commands : List[CommandLineGATK] = Nil
    var dir = ""

    if ( recalOut.getParent != null ) {
      dir = recalOut.getParent+"/"
    }

    var raw_snp = new File(dir+vc.attributes.getProject+".raw_snps.vcf")
    var ug = StandardUnifiedGenotyper(bams, raw_snp)

    commands :+= ug

    var raw_indel = new File(dir+vc.attributes.getProject+".raw_indels.vcf")
    var ig = StandardIndelCalls(bams,raw_indel)

    commands ++= ig

    var prefilt_snp = swapExt(raw_snp,".vcf",".indel_filtered.vcf")
    var iFilt = StandardFilterAtIndels(raw_snp,raw_indel,prefilt_snp)

    commands :+= iFilt

    var annoSNP = prefilt_snp
    if ( refGene != null ) {
      annoSNP = swapExt(prefilt_snp,".vcf",".annotated.vcf")
      var annotate = StandardGenomicAnnotation(prefilt_snp, refGene, annoSNP)

      commands :+= annotate
    }

    var recal = StandardRecalibrateVariants(annoSNP, targetTiTv, recalOut)

    commands ++= recal

    var handFilt = StandardHandfilter(prefilt_snp,handFilteredOut)

    commands :+= handFilt

    return commands
  }

  def StandardCallingPipeline(bams: List[File], indelOut: File, recalOut: File, handFilteredOut: File, targetTiTv: scala.Double, refGene: File = null ) : List[CommandLineGATK] = {
    var commands : List[CommandLineGATK] = Nil
    var dir = ""

    if ( recalOut.getParent != null ) {
      dir = recalOut.getParent+"/"
    }

    var raw_snp = new File(dir+vc.attributes.getProject+".raw_snps.vcf")
    var ug = StandardUnifiedGenotyper(bams, raw_snp)

    commands :+= ug

    var raw_indel = indelOut
    var ig = StandardIndelCalls(bams,raw_indel)

    commands ++= ig

    var prefilt_snp = swapExt(raw_snp,".vcf",".indel_filtered.vcf")
    var iFilt = StandardFilterAtIndels(raw_snp,raw_indel,prefilt_snp)

    commands :+= iFilt

    var annoSNP = prefilt_snp
    if ( refGene != null ) {
      annoSNP = swapExt(prefilt_snp,".vcf",".annotated.vcf")
      var annotate = StandardGenomicAnnotation(prefilt_snp, refGene, annoSNP)

      commands :+= annotate
    }

    var recal = StandardRecalibrateVariants(annoSNP, targetTiTv, recalOut)

    commands ++= recal

    var handFilt = StandardHandfilter(prefilt_snp,handFilteredOut)

    commands :+= handFilt

    return commands
  }
  
  def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

}