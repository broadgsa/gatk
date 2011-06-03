package oneoffs.depristo

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.clipreads.ClippingRepresentation
import scala.io.Source._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction


class PrepareBamsForHomogeneityTesting extends QScript {
  qscript =>
  @Argument(doc="Path to GATK jar",required=false,shortName="gatkjarfile") var gatkJarFile: File = new File("dist/GenomeAnalysisTK.jar")
  @Argument(doc="Path to SamToFastq jar",required=false,shortName="SamToFastqjarfile") var SamToFastqJar: File = new File("/seq/software/picard/current/bin/SamToFastq.jar")
  @Argument(doc="BAM File",required=true, shortName="bam") var bams: File = null
  @Argument(doc="Bundle path",required=false, shortName="bundle") var bundle: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")

//  @Argument(shortName = "R", doc="ref", required=true)
//  var referenceFile: File = bundle + "/human_g1k_v37.fasta"

  @Argument(shortName = "L", doc="intervals", required=false)
  val TARGET_INTERVAL: String = null;

  trait GATKArgs extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    if ( TARGET_INTERVAL != null )
      this.intervalsString = List(TARGET_INTERVAL);
//    this.reference_sequence = referenceFile;
    this.memoryLimit = 2
  }

//  class ClipBAM(@Input in: File, @Output out: File) extends ClipReadsWalker with GATKArgs {
//    this.o = out
//    this.CT = "1-25"
//    this.CR = ClippingRepresentation.HARDCLIP_BASES
//    this.OQ = true
//  }
//
//  case class revert (@Input inBam: File, @Output outBam: File) extends PicardBamFunction {
//    @Output(doc="reverted bam index") var revertedBamIndex = new File(outBam + ".bai")
//    override def inputBams = List(inBAM)
//    override def outputBam = outBam
//    override def commandLine = super.commandLine + " CREATE_INDEX=true "
//    this.isIntermediate = true
//    this.jarFile = qscript.dedupJar
//    this.analysisName = queueLogDir + outBam + ".dedup"
//    this.jobName = queueLogDir + outBam + ".dedup"
//  }

  case class SamToFastq (@Input inBam: File, @Output fastq1: File, @Output fastq2: File, trim: Int) extends JavaCommandLineFunction {
    this.jarFile = qscript.SamToFastqJar
    override def commandLine = super.commandLine +
      Array(
        " INPUT=" + inBam,
        " FASTQ=" + fastq1,
        " VALIDATION_STRINGENCY=SILENT",
        " SECOND_END_FASTQ=" + fastq2,
        " READ1_TRIM=" + trim,
        " READ2_TRIM=" + trim).mkString

    //this.analysisName = queueLogDir + outBam + ".dedup"
    //this.jobName = queueLogDir + outBam + ".dedup"
  }

  def createListFromFile(in: File):List[File] = {
    if (in.toString.endsWith("bam"))
      return List(in)
    var l: List[File] = List()
    for (bam <- fromFile(in).getLines)
      l :+= new File(bam)
    return l
  }

  def script = {
    for ( bam <- createListFromFile(bams) ) {
      val fastq1 = swapExt(bam, ".bam", ".1.fastq")
      val fastq2 = swapExt(bam, ".bam", ".2.fastq")
      add(new SamToFastq(bam, fastq1, fastq2, 25))
    }
  }
}