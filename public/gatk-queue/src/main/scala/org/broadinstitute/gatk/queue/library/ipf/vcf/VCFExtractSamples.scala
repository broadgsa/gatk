/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.library.ipf.vcf

import org.broadinstitute.gatk.queue.function.InProcessFunction
import org.broadinstitute.gatk.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.commandline._
import java.io.{PrintWriter, PrintStream, File}

class VCFExtractSamples(inVCF: File, outVCF: File, samples: List[String]) extends InProcessFunction {
  def this(in: File, out: File, samples: File) = this(in,out, (new XReadLines(samples)).readLines.toList)

  @Input(doc="VCF from which to extract samples") var inputVCF : File = inVCF
  @Output(doc="VCF to which to write the sample-subset vcf") var outputVCF : File = outVCF
  @Argument(doc="The samples to extract from the VCF") var extractSamples : List[String] = samples

  var out : PrintWriter = _
  var columns : List[Int] = (0 to 8).toList

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