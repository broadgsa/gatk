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
import htsjdk.samtools.{SAMSequenceRecord, SAMSequenceDictionary}
import org.broadinstitute.gatk.utils.{GenomeLocParser, GenomeLoc}

class VCFSimpleMerge extends InProcessFunction {
  @Input(doc="VCFs to be merged") var vcfs: List[File] = Nil
  @Input(doc="The reference fasta index") var fai: File = _
  @Output(doc="The final VCF to write to") var outVCF : File = _

  class PeekableXRL(f : File ) {
    var xrl : XReadLines = new XReadLines(f)
    var cur : String = if ( xrl.hasNext ) xrl.next else null

    def hasNext : Boolean = xrl.hasNext || cur != null
    def next : String = {
      var toRet : String = cur
      if ( xrl.hasNext ) {
        cur = xrl.next
      } else {
        cur = null
      }

      return toRet
    }

    def peek : String = cur

  }

  def readHeader( xrl : PeekableXRL ) : List[String] = {
    var toRet : List[String] = Nil
    while ( xrl.hasNext && xrl.peek.startsWith("#") ) {
      toRet :+= xrl.next
    }

    return toRet
  }

  def genomeLoc(xrl : PeekableXRL, p : GenomeLocParser ) : GenomeLoc = {
    var slp = xrl.peek.split("\t",3)
    return p.createGenomeLoc(slp(0),Integer.parseInt(slp(1)),Integer.parseInt(slp(1)))
  }

  def run = {
    var ssd : SAMSequenceDictionary = new SAMSequenceDictionary
    for ( line <- (new XReadLines(fai)).readLines ) {
      val spl = line.split("\\s+")
      val ctig = spl(0)
      val pos = Integer.parseInt(spl(1))
      ssd.addSequence(new SAMSequenceRecord(ctig,pos))
    }

    var xrls : List[PeekableXRL] = vcfs.map( new PeekableXRL(_) )

    var w : PrintWriter = new PrintWriter(new PrintStream(outVCF))

    readHeader(xrls(0)).foreach(u => w.println(u) )

    xrls.foreach(readHeader(_))

    val glp : GenomeLocParser = new GenomeLocParser(ssd)

    var last = ""
    while ( ! xrls.forall( u => ! u.hasNext ) ) {
      val first = xrls.filter( u => u.hasNext).reduceLeft( (a,b) => if ( genomeLoc(a,glp).isBefore(genomeLoc(b,glp))) a else b )
      if ( ! first.peek.equals(last) ) {
        w.println(first.peek)
      }
      last = first.next
    }

    w.close

  }

}