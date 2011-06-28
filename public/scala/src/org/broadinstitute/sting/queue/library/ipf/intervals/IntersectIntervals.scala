package org.broadinstitute.sting.queue.library.ipf.intervals

import org.broadinstitute.sting.queue.function.InProcessFunction
import collection.JavaConversions._
import org.broadinstitute.sting.commandline._
import java.io.{PrintStream, File}
import net.sf.samtools.{SAMSequenceRecord, SAMFileHeader, SAMSequenceDictionary}
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}

class IntersectIntervals(iVals: List[File], outFile: File, bed: Boolean) extends InProcessFunction {
  @Input(doc="List of interval files to find the intersection of") val intervals : List[File] = iVals
  @Output(doc="Output interval file to which to write") val output : File = outFile
  @Argument(doc="Assume the input interval lists are sorted in the proper order") var assumeSorted = false
  @Argument(doc="Is the tdf in bed file (0-based clopen: 0  5 for {1,2,3,4}?") var isBed = bed


  var outStream : PrintStream = _
  var contigs : List[String] = Nil
  var dict : SAMSequenceDictionary = _
  var parser : GenomeLocParser = _

  def run = {
    outStream = new PrintStream(output)
    dict = new SAMSequenceDictionary
    // note: memory hog
    val sources : List[(List[(String,Int,Int)],Int)] = intervals.map(g => asScalaIterator(new XReadLines(g)).map(u => parse(u)).toList).zipWithIndex
    sources.map(u => u._1).flatten.map(u => u._1).distinct.foreach(u => dict.addSequence(new SAMSequenceRecord(u,Integer.MAX_VALUE)))
    parser = new GenomeLocParser(dict)
    sources.map( (u: (List[(String,Int,Int)],Int)) => u._1.map(g => (newGenomeLoc(g),u._2))).flatten.sortWith( (a,b) => (a._1 compareTo b._1) < 0 ).foldLeft[List[List[(GenomeLoc,Int)]]](Nil)( (a,b) => overlapFold(a,b)).map(u => mapIntersect(u)).filter(h => h != null && h.size > 0).foreach(h => writeOut(h))
    outStream.close()
  }

  def writeOut(g : GenomeLoc) : Unit = {
    outStream.print("%s%n".format(g.toString))
  }

  def parse(s : String) : (String,Int,Int) = {
    if ( s.contains(":") ) {
      val split1 = s.split(":")
      val split2 = split1(1).split("-")
      return (split1(0),split2(0).toInt,split2(1).toInt)
    } else {
      val split = s.split("\\s+")
      return (split(0),split(1).toInt + (if(isBed) 1 else 0) ,split(2).toInt - (if(isBed) 1 else 0) )
    }
  }

  def newGenomeLoc(coords : (String,Int,Int) ) : GenomeLoc = {
    parser.createGenomeLoc(coords._1,coords._2,coords._3)
  }

  def overlapFold( a: List[List[(GenomeLoc,Int)]], b: (GenomeLoc,Int) ) : List[List[(GenomeLoc,Int)]] = {
    if ( a.last.forall(u => u._1.overlapsP(b._1)) ) {
      a.init :+ (a.last :+ b)
    } else {
      a :+ ( a.last.dropWhile(u => ! u._1.overlapsP(b._1)) :+ b)
    }
  }

  def mapIntersect( u: List[(GenomeLoc,Int)]) : GenomeLoc = {
    if ( u.map(h => h._2).distinct.sum != range(1,intervals.size).sum ) { // if all sources not accounted for
      null
    }
    u.map(h => h._1).reduceLeft[GenomeLoc]( (a,b) => a.intersect(b) )
  }

  def range(a: Int, b: Int) : Range = new Range(a,b+1,1)

}