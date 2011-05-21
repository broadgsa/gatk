package org.broadinstitute.sting.queue.library.ipf.intervals

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import java.io.{PrintStream, File}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}
import collection.immutable.TreeSet

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

  var intervalCache : TreeSet[GenomeLoc] = _
  val LINES_TO_CACHE : Int = 1000

  def run = {
    output = new PrintStream(outList)
    intervalCache = new TreeSet[GenomeLoc]()(new Ordering[GenomeLoc]{
      def compare(o1: GenomeLoc, o2: GenomeLoc) : Int = { o1.compareTo(o2) }
    })
    parser = new GenomeLocParser(new FastaSequenceFile(ref,true))
    xrl = new XReadLines(inList)
    offsetIn = if (isBed(inType)) 1 else 0
    offsetOut = if( isBed(outType)) 1 else 0
    var line : String = xrl.next
    while ( line.startsWith("@") ) {
      line = xrl.next
    }
    var prevLoc: GenomeLoc = null
    var curLoc: GenomeLoc = null
    var nextLoc : GenomeLoc = parseGenomeInterval(line)
    var linesProcessed : Int = 1
    while ( prevLoc != null || curLoc != null || nextLoc != null ) {
      prevLoc = curLoc
      curLoc = nextLoc
      nextLoc = if ( xrl.hasNext ) parseGenomeInterval(xrl.next) else null
      if ( curLoc != null ) {
        val left: GenomeLoc =  refine(expandLeft(curLoc),prevLoc)
        val right: GenomeLoc =  refine(expandRight(curLoc),nextLoc)
        if ( left != null ) {
          intervalCache += left
        }
        if ( right != null ) {
          intervalCache += right
        }
      }
      linesProcessed += 1
      if ( linesProcessed % LINES_TO_CACHE == 0 ) {
        val toPrint = intervalCache.filter( u => (u.isBefore(prevLoc) && u.distance(prevLoc) > startInt+sizeInt))
        intervalCache = intervalCache -- toPrint
        toPrint.foreach(u => output.print("%s%n".format(repr(u))))
      }
      //System.out.printf("%s".format(if ( curLoc == null ) "null" else repr(curLoc)))
    }

    intervalCache.foreach(u => output.print("%s%n".format(repr(u))))

    output.close()
  }

  def expandLeft(g: GenomeLoc) : GenomeLoc = {
    parser.createGenomeLoc(g.getContig,g.getStart-startInt-sizeInt,g.getStart-startInt)
  }

  def expandRight(g: GenomeLoc) : GenomeLoc = {
    parser.createGenomeLoc(g.getContig,g.getStop+startInt,g.getStop+startInt+sizeInt)
  }

  def refine(newG: GenomeLoc, borderG: GenomeLoc) : GenomeLoc = {
    if ( borderG == null || ! newG.overlapsP(borderG) ) {
      return newG
    } else {
      if ( newG.getStart < borderG.getStart ) {
        if ( borderG.getStart - startInt > newG.getStart ) {
          return parser.createGenomeLoc(newG.getContig,newG.getStart,borderG.getStart-startInt)
        }
      } else {
        if ( borderG.getStop + startInt < newG.getStop ){
          return parser.createGenomeLoc(newG.getContig,borderG.getStop+startInt,newG.getStop)
        }
      }
    }

    null
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
    if ( s.contains(":") ) parser.parseGenomeLoc(s) else parser.createGenomeLoc(sp(0),sp(1).toInt+offsetIn,sp(2).toInt-offsetIn)
  }

  object IntervalFormatType extends Enumeration("INTERVALS","BED","TDF") {
    type IntervalFormatType = Value
    val INTERVALS,BED,TDF = Value

    def convert(s : String) : IntervalFormatType = {
      if ( s.equals("INTERVALS") ) INTERVALS else { if (s.equals("BED") ) BED else TDF}
    }
  }
}