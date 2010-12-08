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
  // TODO -- MAJOR scatter/gather numbers are hard-coded. Really need a function to set it to the right number. These are optimized to get everything on 'hour'

  pm =>

  var stingDirPath : String = stingPath

  class PassFilterAlleles(vcf_files: List[File], out_list: File) extends CommandLineFunction {
    @Input(doc="List of VCFs to extract PF sites from") var vcfs = vcf_files
    @Output(doc="The file to write the site list to") var out_vcf = out_list
    @Argument(doc="Path to the SortByRef script") var sortByRef: String = _
    @Argument(doc="Path to the reference file on disk") var ref: File = _

    def commandLine = {
      "egrep \"FORMAT|format\" %s | cut -f1-8 > %s ; grep PASS %s | tr ':' '\\t' | awk '{print $2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t.\\t.\\t.\"}' | sort -n -k2,2 | uniq | perl %s - %s.fai >> %s".format(
        vcf_files(1).getAbsolutePath, out_vcf.getAbsolutePath, vcf_files.foldLeft[String]("")( (b,a) => b + " " + a.getAbsolutePath), sortByRef, ref.getAbsolutePath, out_vcf.getAbsolutePath
        )
    }
  }

  class Vcf2Intervals(in_vcf: File, out_intervals: File) extends CommandLineFunction {
    @Input(doc="The VCF to convert") var vcf: File = in_vcf
    @Output(doc="The output bed file") var intervals: File = out_intervals

    def commandLine = {
      "grep -v \\\\# %s | awk '{print $1\":\"$2}' | uniq > %s".format(vcf.getAbsolutePath,intervals.getAbsolutePath)
    }
  }

  def MergeBatches( callVCFs: List[File], allBams: List[File], mergedVCF: File, ref: File, size: Int) : List[CommandLineFunction] = {
    var cmds : List[CommandLineFunction] = Nil
    var pfSites : PassFilterAlleles = new PassFilterAlleles(callVCFs,swapExt(mergedVCF,".vcf",".pf.alleles.vcf"))
    pfSites.sortByRef = pm.stingDirPath+"perl/sortByRef.pl"
    pfSites.ref = ref

    cmds :+= pfSites

    var ints : Vcf2Intervals = new Vcf2Intervals(pfSites.out_vcf,swapExt(pfSites.out_vcf,".vcf",".intervals.list"))

    cmds :+= ints
    
    var calcs: List[UGCalcLikelihoods] = allBams.grouped(size).toList.zipWithIndex.map(u => LikelihoodCalc(u._1,ref,ints.intervals,pfSites.out_vcf, new File("batch%d.likelihoods.vcf".format(u._2))))

    cmds ++= calcs

    cmds :+= VariantCallMerge(calcs.map( a => a.out), ref, ints.intervals, mergedVCF)

    return cmds
    
  }
  
  def LikelihoodCalc( bams: List[File], ref: File, intervals: File, alleleVCF: File, outVCF: File ) : UGCalcLikelihoods = {
    var calc = new UGCalcLikelihoods
    calc.input_file ++= bams
    calc.reference_sequence = ref
    calc.jarFile = new File(pm.stingDirPath+"dist/GenomeAnalysisTK.jar")
    calc.downsample_to_coverage = Some(300)
    calc.memoryLimit = if ( bams.size < 5 ) Some(2) else if(bams.size<50) Some(4) else Some(6)
    calc.scatterCount = if (bams.size < 5  ) 1 else if (bams.size < 50) 60 else 120
    calc.min_base_quality_score = Some(22)
    calc.min_mapping_quality_score = Some(20)
    calc.genotype = true
    calc.output_all_callable_bases = true
    calc.out = outVCF
    calc.rodBind :+= new RodBind("allele","VCF",alleleVCF)
    calc.intervals :+= intervals

    return calc
  }

  def VariantCallMerge( likelihoods: List[File], ref: File, intervals: File, output: File) : UGCallVariants = {
    var call = new UGCallVariants
    call.reference_sequence = ref
    call.jarFile = new File(pm.stingDirPath+"dist/GenomeAnalysisTK.jar")
    call.intervals :+= intervals
    call.memoryLimit = Some(8)
    call.out = output
    call.rodBind ++= likelihoods.map( a => new RodBind("variant"+a.getName.replace(".vcf",""),"vcf",a) )
    call.scatterCount = 30

    return call
  }

  def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

}