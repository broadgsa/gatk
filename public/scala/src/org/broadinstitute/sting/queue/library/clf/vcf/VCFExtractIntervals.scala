package org.broadinstitute.sting.queue.library.clf.vcf

import java.io.File
import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import org.broadinstitute.sting.queue.function.CommandLineFunction

class VCFExtractIntervals(inVCF: File, outList: File, passOnly: Boolean ) extends CommandLineFunction {
  def this(vcf: File, list: File ) = this(vcf,list,false)
  def this(vcf: File) = this(vcf, new File(vcf.getAbsolutePath.replace("vcf","intervals.list")),false)
  def this(vcf: File, removeFilters: Boolean) = this(vcf, new File(vcf.getAbsolutePath.replace("vcf","intervals.list")), removeFilters)

  @Input(doc="The VCF from which to extract an interval list") var inputVCF : File = inVCF
  @Output(doc="The file to write the interval list to") var outputList : File = outList
  @Argument(doc="Whether to use all sites, or only the unfiltered sites") var usePFOnly: Boolean = passOnly

  def commandLine = {
    if ( usePFOnly ) "grep PASS %s | awk '{print $1\":\"$2}' | uniq > %s".format(inputVCF.getAbsolutePath,outputList.getAbsolutePath)
    else "grep -v \\\\# %s | awk '{print $1\":\"$2}' | uniq > %s".format(inputVCF.getAbsolutePath,outputList.getAbsolutePath)
  }

}