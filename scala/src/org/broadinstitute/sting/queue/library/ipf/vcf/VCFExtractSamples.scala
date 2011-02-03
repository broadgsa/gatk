package org.broadinstitute.sting.queue.library.ipf.vcf

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.sting.commandline._
import java.io.{PrintWriter, PrintStream, File}

class   VCFExtractSamples(inVCF: File, outVCF: File, samples: List[String]) extends InProcessFunction {
  def this(in: File, out: File, samples: File) = this(in,out, (new XReadLines(samples)).readLines.toList)

  @Input(doc="VCF from which to extract samples") var inputVCF : File = inVCF
  @Output(doc="VCF to which to write the sample-subset vcf") var outputVCF : File = outVCF
  @Argument(doc="The samples to extract from the VCF") var extractSamples : List[String] = samples

  var out : PrintWriter = _
  var columns : List[Int] = 0 to 8 toList

  def run = {
    out = new PrintWriter(new PrintStream(outputVCF))
    asScalaIterator(new XReadLines(inputVCF)).foreach(subset)
    out.close
  }

  def subset( line : String ) {
    if ( line.startsWith("##") ) {
      out.print("%s%n".format(line))
    } else {
      val spline = line.split("\t")
      if ( spline(0).equals("#CHROM") ) {
        columns ++= spline.zipWithIndex.filter( p => samples.contains(p._1) ).map( p => p._2 )
      }

      out.print("%s%n".format(columns.map(p => spline(p)).reduceLeft(_ + "\t" + _)))
    }
  }
}