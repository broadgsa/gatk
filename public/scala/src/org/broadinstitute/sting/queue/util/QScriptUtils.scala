package org.broadinstitute.sting.queue.util

import java.io.File
import io.Source._
import net.sf.samtools.{SAMReadGroupRecord, SAMFileReader}

import collection.JavaConversions._


/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/14/11
 * Time: 4:57 PM
 * To change this template use File | Settings | File Templates.
 */

object QScriptUtils {

  /**
   * Takes a bam list file and produces a scala list with each file allowing the bam list
   * to have empty lines and comment lines (lines starting with #).
   */
  def createListFromFile(in: File):List[File] = {
    // If the file provided ends with .bam, .fasta or .fq, it is not a bam list, we treat it as a single file.
    // and return a list with only this file.
    if (in.toString.endsWith(".bam") || in.toString.endsWith(".fasta") || in.toString.endsWith(".fq"))
      return List(in)

    var list: List[File] = List()
    for (file <- fromFile(in).getLines)
      if (!file.startsWith("#") && !file.isEmpty )
        list :+= new File(file.trim())
    list.sortWith(_.compareTo(_) < 0)
  }

  /**
   * Returns the number of contigs in the BAM file header.
   */
  def getNumberOfContigs(bamFile: File): Int = {
    val samReader = new SAMFileReader(bamFile)
    samReader.getFileHeader.getSequenceDictionary.getSequences.size()
  }

  /**
   * Check if there are multiple samples in a BAM file
   */
  def hasMultipleSamples(readGroups: java.util.List[SAMReadGroupRecord]): Boolean = {
    var sample: String = ""
    for (r <- readGroups) {
      if (sample.isEmpty)
        sample = r.getSample
      else if (sample != r.getSample)
          return true;
    }
    false
  }


  def ?[A <: AnyRef](ref: A): Option[A] =
    if (ref eq null) None else Some(ref)
}