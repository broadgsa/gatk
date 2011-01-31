package org.broadinstitute.sting.queue.library.ipf.vcf

import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{PrintStream, PrintWriter, File}


class VCFExtractSites( vcf: File, output: File) extends InProcessFunction {
  @Input(doc="VCF file from which to extract sites") var inVCF: File = vcf
  @Output(doc="Sites VCF file to write to") var outVCF: File = output
  @Argument(doc="Keep non-PASS sites") var keepFilters: Boolean = false
  @Argument(doc="Keep info field") var keepInfo : Boolean = true
  @Argument(doc="Keep qual field") var keepQual : Boolean = true

  def lineMap( line: String ) : String = {
    if ( line.startsWith("##") ) { return line }
    val spline = line.split("\t",9)

    if ( spline(6) == "PASS" || keepFilters ) {
      if ( ! keepInfo ) {
        spline(7) = "."
      }
      if ( ! keepQual ) {
        spline(5) = "."
      }
      return spline.slice(0,8).reduceLeft( _ + "\t" + _ )
    }

    return ""
  }

  def run {
    var w: PrintWriter = new PrintWriter( new PrintStream(outVCF) )
    ( new XReadLines(inVCF) ).readLines().map(lineMap).view.filter( u => u != "" ).foreach( u => w.println(u) )
  }

}