package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util._
import java.io.File
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.queue.function.CommandLineFunction

class ProjectManagement(stingPath: String) {

  pm =>

  var stingDirPath : String = stingPath

  class PassFilterSites(vcf_files: List[File], out_list: File) extends CommandLineFunction {
    @Input(doc="List of VCFs to extract PF sites from") var vcfs = vcf_files
    @Output(doc="The file to write the site list to") var out_intervals = out_list
    @Argument(doc="Path to the SortByRef script") var sortByRef: String = _
    @Argument(doc="Path to the reference file on disk") var ref: File = _

    def commandLine = {
      "grep PASS %s | tr ':' '\\t' | awk '{print $2\"\\t\"$3}' | sort -n -k2,2 | uniq | perl %s - %s.fai | awk '{print $1\":\"$2}' > %s".format(
        vcf_files.foldLeft[String]("")( (b,a) => b + " " + a.getAbsolutePath), sortByRef, ref.getAbsolutePath, out_list.getAbsolutePath
        )
    }
  }

  def MergeBatches( callVCFs: List[File], allBams: List[File], mergedVCF: File, ref: File) : List[CommandLineFunction] = {
    var cmds : List[CommandLineFunction] = Nil
    var pfSites : PassFilterSites = new PassFilterSites(callVCFs,swapExt(mergedVCF,".vcf",".pf.intervals.list"))
    pfSites.sortByRef = pm.stingDirPath+"perl/sortByRef.pl"
    pfSites.ref = ref

    cmds :+= pfSites
    
    var calcs: List[UGCalcLikelihoods] = allBams.map( a => LikelihoodCalc(a,ref,pfSites.out_intervals) )

    cmds ++= calcs

    cmds :+= VariantCallMerge(calcs.map( a => a.out), ref, pfSites.out_intervals, mergedVCF)

    return cmds
    
  }

  def LikelihoodCalc( bam: File, ref: File, intervals: File ) : UGCalcLikelihoods = {
    var calc = new UGCalcLikelihoods
    calc.input_file :+= bam
    calc.reference_sequence = ref
    calc.jarFile = new File(pm.stingDirPath+"dist/GenomeAnalysisTK.jar")
    calc.intervals :+= intervals
    calc.downsample_to_coverage = Some(300)
    calc.memoryLimit = Some(2)
    calc.min_base_quality_score = Some(22)
    calc.min_mapping_quality_score = Some(20)
    calc.genotype = true
    calc.output_all_callable_bases = true
    calc.out = swapExt(bam,".bam",".likelihoods.vcf")

    return calc
  }

  def VariantCallMerge( likelihoods: List[File], ref: File, intervals: File, output: File) : UGCallVariants = {
    var call = new UGCallVariants
    call.reference_sequence = ref
    call.jarFile = new File(pm.stingDirPath+"dist/GenomeAnalysisTK.jar")
    call.intervals :+= intervals
    call.memoryLimit = Some(4)
    call.out = output
    call.rodBind ++= likelihoods.map( a => new RodBind("variant"+a.getName.replace(".vcf",""),"vcf",a) )

    return call
  }

  def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

}