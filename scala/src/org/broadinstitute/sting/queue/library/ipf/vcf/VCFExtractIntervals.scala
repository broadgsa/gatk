package org.broadinstitute.sting.queue.library.ipf.vcf

import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{PrintStream, PrintWriter, File}

class VCFExtractIntervals(inVCF: File, outList: File, useFilterSites: Boolean) extends InProcessFunction {
  def this(in : File, out: File) = this(in,out,true)
  def this(in : File) = this(in,new File(in.getAbsolutePath.replace(".vcf",".intervals.list")),true)

  @Input(doc="The VCF to convert to an interval list") var vcfIn : File = inVCF
  @Output(doc="The intervals file to write to") var listOut : File = outList
  @Argument(doc="Keep filtered sites?") var keepFilters : Boolean = useFilterSites

  var out : PrintWriter = _

  def run = {
    out = new PrintWriter(new PrintStream(listOut))
    asScalaIterator(new XReadLines(vcfIn)).foreach(vcf2int)
  }

  def vcf2int( vcfLine: String ) : Unit = {
    var spline = vcfLine.split("\t")
    if ( ! vcfLine.startsWith("#") && (spline(6).equals("PASS") || keepFilters) ) {
      out.print("%s:%s%n".format(spline(0),spline(1)))
    }
  }

}