package org.broadinstitute.sting.queue.library.ipf

import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{PrintStream, PrintWriter, File}
import collection.immutable.HashMap


class SortByRef( input: File, reference: File, output: File ) extends InProcessFunction {
  @Input(doc="The file to be sorted") var inFile: File = input
  @Input(doc="The reference fasta index") var fai: File = reference
  @Output(doc="The file to write the sorted file to") var outFile :  File = output
  @Argument(doc="The character or expression that separates entries") var separator : String = "\t"
  @Argument(doc="The position of the contig in the file (1-based)") var pos: Int = 1
  @Argument(doc="Comment characters (lines will be ignored)") var comment: List[String] = List("#")

  var contigMap: HashMap[String,Int] = new HashMap[String,Int];

  def contigVal( line : String ) : Int = {
    if ( comment.contains(line.charAt(0)) ) {
      return -1;
    }

    if ( contigMap.size < 1 ) { // no contigs
      ( new XReadLines(fai)).readLines.map( u => u.split("\t").head).zipWithIndex.foreach( u => contigMap += u )
    }

    return contigMap( line.split(separator)(pos-1) )
  }

  def run = {
    var w : PrintWriter = new PrintWriter(new PrintStream(outFile))
    ( new XReadLines(inFile) ).readLines.sortBy(contigVal).foreach( u => w.println(u) )    
  }
}