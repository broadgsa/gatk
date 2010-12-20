package org.broadinstitute.sting.queue.library.ipf.intervals

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{PrintStream, File}

/**
 * Utility for explicitly grouping intervals together in sets of size num
 * @Author chartl
 */

class GroupIntervals(iList: File, num: Int, head: Boolean, dir: Option[String]) extends InProcessFunction {
  @Input(doc="Interval list to split into groups") var intervalList : File = iList
  @Output(doc="The output interval lists created") var outputList : List[File] = Nil
  @Argument(doc="The number of groups wanted") var nGroups : Int = num
  @Argument(doc="The base location you want the files") var directory : Option[String] = dir
  @Argument(doc="Has header (@ lines)") var hasHeader : Boolean = head

  var written : Int = 0
  var header : List[String] = Nil

  def run = {
    written = 0
    header = Nil
    asScalaIterator(new XReadLines(intervalList)).filter(u => ! isHeader(u)).grouped(num).foreach(n => write(n))
  }

  def write(lines : Seq[String]) = {
    val oFile : File = if ( dir.isEmpty) new File(intervalList.getAbsolutePath+".group%d".format(written))else new File(dir.get+"/"+intervalList.getAbsolutePath+".group%d".format(written))
    outputList :+= oFile
    val output : PrintStream = new PrintStream(oFile)
    header.foreach( u => output.print("%s%n".format(u)))
    lines.foreach(u => output.print("%s%n".format(u)))
    output.close()
  }

  def isHeader(s : String) : Boolean = {
    if ( hasHeader && s.startsWith("@") ) {
      header :+= s
      return true
    } else {
      return false
    }
  }

}