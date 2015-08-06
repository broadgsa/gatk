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

import collection.JavaConversions._
import org.broadinstitute.gatk.queue.function.InProcessFunction
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.utils.text.XReadLines
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
    var elems = asScalaIterator(new XReadLines(vcfIn)).map(vcf2int).filter(p => !p.equals(""))
    var prev : String = null
    if ( elems.hasNext ) {
      prev = elems.next
    }
    var cur : String = null
    if ( elems.hasNext ) {
      cur = elems.next
      while ( elems.hasNext ) {
        out.printf("%s%n",prev)
        while ( cur.equals(prev) && elems.hasNext && !cur.equals("") ) {
          cur = elems.next
        }

        if ( ! cur.equals(prev) ) {
          if ( elems.hasNext ) {
            prev = cur
            cur = elems.next
          }
        }
      }
      out.printf("%s%n",prev)
      out.printf("%s%n",cur)
    } else {
      out.printf("%s%n",prev)
    }

    out.close
  }

  def vcf2int( vcfLine: String ) : String = {
    var spline = vcfLine.split("\t")
    if ( ! vcfLine.startsWith("#") && (spline(6).equals("PASS") || spline(6).equals(".") || keepFilters) ) {
      return("%s:%s".format(spline(0),spline(1)))
    } else {
      return ""
    }
  }

}