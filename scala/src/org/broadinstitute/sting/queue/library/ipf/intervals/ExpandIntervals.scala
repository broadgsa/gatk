package org.broadinstitute.sting.queue.library.ipf.intervals

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import java.io.{PrintStream, File}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}

// todo -- this is unsafe. Need to use a reference dictionary to ensure no off-contig targets are created
class ExpandIntervals(in : File, start: Int, size: Int, out: File, ref: File, ipType: String, opType: String) extends InProcessFunction {
  @Input(doc="The interval list to expand") val inList : File = in
  @Input(doc="The reference sequence") val refDict : File = ref
  @Argument(doc="Number of basepair to start the expanded interval") val startInt : Int = start
  @Argument(doc="Number of baispair to stop the expanded interval") val sizeInt : Int = size
  @Output(doc="The output intervals file to write to") val outList : File = out
  @Argument(doc="The output format for the intervals") val outTypeStr = opType
  @Argument(doc="The input format for the intervals") val inTypeStr = ipType

  var output : PrintStream = _
  var parser : GenomeLocParser = _
  var xrl : XReadLines = _
  val outType = IntervalFormatType.convert(outTypeStr)
  val inType = IntervalFormatType.convert(inTypeStr)

  var offsetIn : Int = 0
  var offsetOut : Int = 0

  var first : Boolean = true
  var lastTwo : (GenomeLoc,GenomeLoc) = _

  def run = {
    first = true
    lastTwo = null
    output = new PrintStream(outList)
    parser = new GenomeLocParser(new FastaSequenceFile(ref,true))
    xrl = new XReadLines(inList)
    offsetIn = if (isBed(inType)) 1 else 0
    offsetOut = if (isBed(outType)) 1 else 0
    asScalaIterable(xrl).filter( ! _.startsWith("@")).map(parseGenomeInterval(_)).sliding(3).map(a => a.toList).foreach(a => expand(a(0),a(1),a(2)))
    expand(lastTwo._1,lastTwo._2,null)
    output.close()
  }

  def expand(previous: GenomeLoc, current: GenomeLoc, next: GenomeLoc) : Unit = {
    if ( first ) {
      first = false
      expand(null,previous,current)
    }
    lastTwo = (current,next)

    if ( start > 0 ) {
      val new1 = parser.createGenomeLoc(current.getContig,current.getStart-startInt-sizeInt,current.getStart-startInt)
      val new2 = parser.createGenomeLoc(current.getContig,current.getStop+startInt,current.getStop+startInt+sizeInt)
      if ( ok(new1,previous,next) ) {
        //System.out.println("Printing! %s".format(repr(new1)))
        output.print("%s%n".format(repr(new1)))
      }
      if ( ok(new2,previous,next) ) {
        //System.out.println("Printing! %s".format(repr(new2)))
        output.print("%s%n".format(repr(new2)))
      }
    } else {
      output.print("%s%n".format(repr(parser.createGenomeLoc(current.getContig,current.getStart-sizeInt,current.getStop+sizeInt))))
    }
  }

  def ok( loc : GenomeLoc, prev: GenomeLoc, next: GenomeLoc  ) : Boolean = {
    //System.out.println("%s - %s - %s".format(repr(next),repr(loc),repr(previous)))
    ( next == null || loc.distance(next) >= start) && (prev == null || loc.distance(prev) >= start)
  }

  def repr(loc : GenomeLoc) : String = {
    if ( loc == null ) return "null"
    if ( outType == IntervalFormatType.INTERVALS ) {
      return "%s:%d-%d".format(loc.getContig,loc.getStart,loc.getStop)
    } else {
      return "%s\t%d\t%d".format(loc.getContig,loc.getStart-offsetOut,loc.getStop+offsetOut)
    }
  }

  def isBed(t: IntervalFormatType.IntervalFormatType) : Boolean = {
    t == IntervalFormatType.BED
   }

  def parseGenomeInterval( s : String ) : GenomeLoc = {
    val sp = s.split("\\s+")
    // todo -- maybe specify whether the bed format [0,6) --> (1,2,3,4,5) is what's wanted  
    if ( s.contains(":") ) parser.parseGenomeInterval(s) else parser.createGenomeLoc(sp(0),sp(1).toInt+offsetIn,sp(2).toInt-offsetIn)
  }

  object IntervalFormatType extends Enumeration("INTERVALS","BED","TDF") {
    type IntervalFormatType = Value
    val INTERVALS,BED,TDF = Value

    def convert(s : String) : IntervalFormatType = {
      if ( s.equals("INTERVALS") ) INTERVALS else { if (s.equals("BED") ) BED else TDF}
    }
  }
}