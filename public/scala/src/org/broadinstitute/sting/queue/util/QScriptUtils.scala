/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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
 */

object QScriptUtils {

  /**
   * Takes a bam list file and produces a scala sequence with each file allowing the bam list
   * to have empty lines and comment lines (lines starting with #).
   */
  def createSeqFromFile(in: File):Seq[File] = {
    // If the file provided ends with .bam, .fasta, fastq or .fq, it is not a bam list, we treat it as a single file.
    // and return a list with only this file.
    if (in.toString.toUpperCase.endsWith(".BAM")   ||
        in.toString.toUpperCase.endsWith(".FASTA") ||
        in.toString.toUpperCase.endsWith(".FQ")    ||
        in.toString.toUpperCase.endsWith("FASTQ") )
      return Seq(in)

    var list: Seq[File] = Seq()
    for (file <- fromFile(in).getLines())
      if (!file.startsWith("#") && !file.isEmpty )
        list :+= new File(file.trim())
//    list.sortWith(_.compareTo(_) < 0)
    list
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
}
