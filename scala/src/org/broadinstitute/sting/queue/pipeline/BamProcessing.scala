package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.queue.extensions.gatk._
import java.io.File
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import net.sf.picard.reference.ReferenceSequenceFileFactory
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}
import org.broadinstitute.sting.utils.interval.IntervalUtils
import collection.mutable.{ListBuffer, HashMap}
import collection.JavaConversions
import java.util.Arrays
import org.broadinstitute.sting.queue.util.{PipelineUtils, IOUtils}
import org.broadinstitute.sting.commandline.{Output, Input}
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction

class BamProcessing(yaml: File, gatkJar: File, fixMatesJar: File) {
  library =>

  var attributes: Pipeline = YamlUtils.load(classOf[Pipeline],yaml)

  trait StandardCommandLineGATK extends CommandLineGATK {
    this.reference_sequence = library.attributes.getProject.getReferenceFile
    this.intervals = library.attributes.getProject.getIntervalList
    this.DBSNP = library.attributes.getProject.getDbsnpFile
    this.memoryLimit = Some(2)
    this.jarFile = library.gatkJar
  }

  /**
   * @Doc: Creates a standard realigner target creator CLF given a bam, an output file, and the contigs over which to run
   * @Returns: A CLF for the realigner target creator job
   */
  protected def StandardRealignerTargetCreator(bam: File, contigs: List[String], output: File) : RealignerTargetCreator = {
    var rtc = new RealignerTargetCreator with StandardCommandLineGATK
    rtc.intervals = null
    rtc.intervalsString = contigs
    rtc.input_file :+= bam
    rtc.out = output
    rtc.analysisName = "RealignerTargetCreator"

    return rtc
  }

  /**
   * @Doc: Creates a standard indel cleaner CLF given a bam, the results of the target creator, and an output .bam file
   * @Returns: A CLF for the indel cleaning job
   */
  protected def StandardIndelCleaner(bam: File, contigs: List[String], targets: File, outBam: File) : IndelRealigner = {
    var realigner = new IndelRealigner with StandardCommandLineGATK
    realigner.intervalsString = contigs
    realigner.intervals = null
    realigner.input_file :+= bam
    realigner.out = outBam
    realigner.targetIntervals = targets
    realigner.analysisName = "IndelClean"
    realigner.bam_compression = Some(0)

    return realigner    
  }

  /**
   * @Doc: Creates a standard split-by-contig indel cleaner job for a given bam file, RTC output, and bam to merge everything to
   * @Returns: A list of CLFs (todo -- wrapped in a Pipeline)
   */
  protected def StandardIndelCleanBam(bam: File, jobContigs: List[List[String]], targets: File, cleanedBam: File) : List[CommandLineFunction] = {
    var cmds : List[CommandLineFunction] = Nil
    var jobSpecs : List[(File,File,List[String])] = jobContigs.map[(File,File,List[String]),List[(File,File,List[String])]](
      ctigs => { (bam, swapExt(bam,".bam",".%s.bam".format(ctigs.mkString("_"))), ctigs) }
      )
    var bamsToMerge : List[File] = Nil
    for ( spec <- jobSpecs ) {
      cmds :+= StandardIndelCleaner(spec._1,spec._3,targets,spec._2)
      bamsToMerge :+= spec._2
    }

    cmds :+= StandardPicardFixMates(bamsToMerge,cleanedBam,library.fixMatesJar)

    return cmds

  }

  /**
   * @Doc: Given a list of (pairs of) bams and cleaned bams to write to, and a number of jobs, creates a set of
   * command line functions to do the target-creating, splitting, cleaning, and merging, returning that list
   * of command line functions
   * @Returns: A list of command line functions for the full indel realignment pipeline from the collection
   * of uncleaned bams to the collection of cleaned bams
   */
  protected def StandardIndelRealign( bamsUncleanCleanPairs: List[(File,File)], nJobs: Int = 1 ) : List[CommandLineFunction] = {
    val contigsForJobs : List[List[String]] = PipelineUtils.smartSplitContigs(library.attributes.getProject.getReferenceFile, library.attributes.getProject.getIntervalList, nJobs)
    var commands : List[CommandLineFunction] = Nil
    for ( bamPair <- bamsUncleanCleanPairs ) {
      val rtc : RealignerTargetCreator = StandardRealignerTargetCreator(bamPair._1,contigsForJobs.foldLeft[List[String]](Nil)( (a,b) => a ::: b), swapExt(bamPair._1,".bam",".targets") )
      val icbs : List[CommandLineFunction] = StandardIndelCleanBam(bamPair._1,contigsForJobs,rtc.out,bamPair._2)
      val sam : SamtoolsIndexFunction = new SamtoolsIndexFunction
      sam.bamFile = bamPair._2
      sam.analysisName = "SamtoolsIndex"
      commands :+= rtc
      commands ++= icbs
      commands :+= sam
    }

    return commands
  }

  /**
   * @Doc: Merges N bam files into one bam file, fixing mate pairs in the process; does not assume they are sorted
   * @Returns: Command line function for the merge, fix-mate, and sort operation
   */
  protected def StandardPicardFixMates(inBams: List[File], outBam: File, picardJar: File) : CommandLineFunction = {
    var pfm : PicardFixMates = new PicardFixMates
    pfm.bams = inBams
    pfm.outBam = outBam
    pfm.jarFile = picardJar
    pfm.assumeSorted = Some(false)
    pfm.memoryLimit = Some(4)
    pfm.analysisName = "FixMates"

    return pfm
  }

  class PicardFixMates extends PicardBamJarFunction {
    @Input(doc="input bam files") var bams: List[File] = Nil
    @Output(doc="output bam file") var outBam: File = null

    def inputBams: List[File] = bams
    def outputBam: File = outBam

  }


  def swapExt(file: File, oldExtension: String, newExtension: String) =
     new File(file.getName.stripSuffix(oldExtension) + newExtension)
  
}
