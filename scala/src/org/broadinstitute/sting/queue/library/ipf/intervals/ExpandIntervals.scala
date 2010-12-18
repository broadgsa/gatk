package org.broadinstitute.sting.queue.library.ipf

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import java.io.{PrintStream, File}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}

// todo -- this is unsafe. Need to use a reference dictionary to ensure no off-contig targets are created
class ExpandIntervals(in : File, start: Int, size: Int, out: File, ref: File, opType: String) extends InProcessFunction {
  @Input(doc="The interval list to expand") val inList : File = in
  @Input(doc="The reference sequence") val refDict : File = ref
  @Argument(doc="Number of basepair to start the expanded interval") val startInt : Int = start
  @Argument(doc="Number of baispair to stop the expanded interval") val sizeInt : Int = size
  @Output(doc="The output intervals file to write to") val outList : File = out
  @Argument(doc="The output format for the intervals") val outTypeStr = opType

  val output : PrintStream = new PrintStream(outList)
  val parser : GenomeLocParser = new GenomeLocParser(new FastaSequenceFile(ref,true))
  val xrl : XReadLines = new XReadLines(inList)
  val outType = if ( outTypeStr.equals("INTERVALS") ) IntervalOutputType.INTERVALS else IntervalOutputType.BED

  var previous : GenomeLoc = _
  var current : GenomeLoc = _
  var next : GenomeLoc = _

  def run = {
    // todo -- there's got to be a special view that gives prevous, current, next
    asScalaIterable(xrl).filter( ! _.startsWith("@")).map(parseGenomeInterval(_)).foreach(expand(_))
  }

  def expand(loc : GenomeLoc) : Unit = {
    if ( start > 0 ) {
      previous = current
      current = next
      next = loc
      if ( current == null ) {
        return
      }
      val new1 = parser.createGenomeLoc(current.getContig,current.getStart-startInt-sizeInt,current.getStart-startInt)
      val new2 = parser.createGenomeLoc(current.getContig,current.getStop+startInt,current.getStop+startInt+sizeInt)
      if ( ok(new1) ) {
        //System.out.println("Printing! %s".format(repr(new1)))
        output.print("%s%n".format(repr(new1)))
      }
      if ( ok(new2) ) {
        //System.out.println("Printing! %s".format(repr(new2)))
        output.print("%s%n".format(repr(new2)))
      }
      previous = current
    } else {
      output.print("%s%n".format(repr(parser.createGenomeLoc(loc.getContig,loc.getStart-sizeInt,loc.getStop+sizeInt))))
    }
  }

  def ok( loc : GenomeLoc ) : Boolean = {
    //System.out.println("%s - %s - %s".format(repr(next),repr(loc),repr(previous)))
    ( next == null || loc.distance(next) >= start) && (previous == null || loc.distance(previous) >= start)
  }

  def repr(loc : GenomeLoc) : String = {
    if ( loc == null ) return "null"
    if ( outType == IntervalOutputType.INTERVALS ) {
      return "%s:%d-%d".format(loc.getContig,loc.getStart,loc.getStop)
    } else if ( outType == IntervalOutputType.BED ) {
      return "%s\t%d\t%d".format(loc.getContig,loc.getStart-1,loc.getStop)
    } else {
      return "FORMAT?"
    }
  }

  def parseGenomeInterval( s : String ) : GenomeLoc = {
    val sp = s.split("\\s+")
    if ( s.contains(":") ) parser.parseGenomeInterval(s) else parser.createGenomeLoc(sp(0),sp(1).toInt,sp(2).toInt)
  }

  object IntervalOutputType extends Enumeration("INTERVALS","BED") {
    type IntervalOutputType = Value
    val INTERVALS,BED = Value
  }
}