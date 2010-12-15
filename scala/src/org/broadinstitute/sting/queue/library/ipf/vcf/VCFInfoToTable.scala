package org.broadinstitute.sting.queue.library.ipf.vcf

import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._
import java.io.{PrintStream, PrintWriter, File}
import collection.immutable.HashSet
import collection.mutable.HashMap
import org.broadinstitute.sting.queue.function.InProcessFunction

class VCFInfoToTable(vcf: File, table: File, annots: List[String], keepFilter: Boolean) extends InProcessFunction {
  def this(in: File, out: File, anns: List[String]) = this(in,out,anns,true)
  def this(in: File, out: File) = this(in,out,Nil,true)
  def this(in: File) = this(in, new File(in.getAbsolutePath.replace(".vcf",".info.table")), Nil, true)
  def this(in: File, anns: List[String]) = this(in, new File(in.getAbsolutePath.replace(".vcf",".info.table")),anns,true)
  def this(in: File, anns: java.util.List[String]) = this(in, new File(in.getAbsolutePath.replace(".vcf",".info.table")),anns.toList,true)
  def this(in: File, anns: List[String], keep: Boolean) = this(in, new File(in.getAbsolutePath.replace(".vcf",".info.table")),anns,keep)
  def this(in: File, anns: java.util.List[String], keep: Boolean) = this(in, new File(in.getAbsolutePath.replace(".vcf",".info.table")),anns.toList,keep)
  def this(in: File, out: File, anns: java.util.List[String], keep: Boolean) = this(in,out,anns.toList,keep)

  @Input(doc="VCF file from which to extract annotion") var inVCF: File = vcf
  @Output(doc="Table file to which to write") var outTable: File = table
  @Argument(doc="Annotations to extract from info field") var annotations: List[String] = annots
  @Argument(doc="Keep filtered records?") var keepFilteredRecs: Boolean = keepFilter

  // set as vars so pipelnes can hack these values
  var NO_KEY : String = "NA"
  var PF_KEY : String = "PASS"

  var out : PrintWriter = _
  var annotation_set : HashSet[String] = new HashSet[String]

  //todo -- Khalid: Why is run not being called?
  def run = {
    logger.debug("RUN IS CALLED")
    annotation_set ++= annotations
    out = new PrintWriter(new PrintStream(outTable))
    asScalaIterator(new XReadLines(inVCF)).foreach(lineToTable)
  }

  def lineToTable(line : String) = {
    if ( ! line.startsWith("#") ) {
      val spline = line.split("\t")
      if ( spline(6).equals(PF_KEY) || keepFilteredRecs ) {
        val iMap = spline(7).split(";").map(_.split("=")).filter(p => annotation_set.contains(p.apply(0))).foldLeft(new HashMap[String,String])( (a,b) => a += new Tuple2(b.apply(0),b.apply(1)) )
        out.print("%s%n".format(annotations.map( u => {
          if ( iMap.contains(u) ) {
            iMap.get(u)
          } else {
            NO_KEY
          }
        }).reduceLeft((a,b) => a + "\t" + b)))
      }
    }
  }
}