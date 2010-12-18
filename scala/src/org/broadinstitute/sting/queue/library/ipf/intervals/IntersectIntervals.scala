package org.broadinstitute.sting.queue.library.ipf.intervals

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import java.io.{PrintStream, File}
import net.sf.samtools.{SAMSequenceRecord, SAMFileHeader, SAMSequenceDictionary}
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}

abstract class IntersectIntervals(iVals: List[File], outFile: File) extends InProcessFunction {
 /* @Input(doc="List of interval files to find the intersection of") val intervals : List[File] = iVals
  @Output(doc="Output interval file to which to write") val output : File = outFile
  @Argument(doc="List of contigs for the interval list (defaults to human b37 autosome); IN ORDER")
  var contigs : List[String] = (new Range(1,22)).map(u => "%d".format(u)) ::: List("X","Y","MT")
  @Argument(doc="Assume the input interval lists are sorted in the proper order") var assumeSorted = false


  val out : PrintStream = new PrintStream(output)

  def run = {
    var dict : SAMSequenceDictionary = new SAMSequenceDictionary
    contigs.map( u => new SAMSequenceRecord(u, Integer.MAX_VALUE ) ).foreach(u => dict.addSequence(u))
    val parser : GenomeLocParser = new GenomeLocParser(dict)
    val locList : IntervalIterator = new IntervalIterator(intervals,parser,assumeSorted)
    
  }

  class IntervalIterator(ivals : List[File], parser: GenomeLocParser, sorted: Boolean) extends Iterable[GenomeLoc] {
    val readers : List[XReadLinesBuffer] = ivals.map(f => new XReadLinesBuffer(new XReadLines(f), parser))
    val filesAreSorted : Boolean = sorted

    def hasNext = ! readers.filter( _.hasNext ).isEmpty
    def next = if ( assumeSorted ) nextElement else nextSortedElement
    def nextElement =


  }

  class XReadLinesBuffer( x: XReadLines, p: GenomeLocParser) extends Iterable[GenomeLoc] {
    val xrl: XReadLines = x
    val parser : GenomeLocParser = p
    var current : GenomeLoc = _

    def hasNext : Boolean = curLine != null || x.hasNext
    def next : String = {
      if ( curLine == null && x.hasNext ) curLine = x.next
      val retLine = curLine
      if ( x.hasNext ) curLine = x.next else curLine = null
      retLine
    }

  } */

}