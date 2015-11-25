/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.util

import java.io.File
import org.apache.commons.io.FilenameUtils
import scala.io.Source._
import htsjdk.samtools.SAMFileReader
import htsjdk.variant.vcf.{VCFHeader, VCFCodec}
import scala.collection.JavaConversions._
import htsjdk.tribble.AbstractFeatureReader

object VCF_BAM_utilities {

  def getSamplesFromVCF(vcfFile: File): List[String] = {
    AbstractFeatureReader.getFeatureReader(vcfFile.getPath, new VCFCodec()).getHeader.asInstanceOf[VCFHeader].getGenotypeSamples.toList
  }

  def getSamplesInBAM(bam: File): List[String] = {
    return new SAMFileReader(bam).getFileHeader().getReadGroups().toList.map(srgr => srgr.getSample()).toSet.toList
  }

  def parseBAMsInput(bamsIn: File): List[File] = FilenameUtils.getExtension(bamsIn.getPath) match {
    case "bam" => return List(bamsIn)
    case "list" => return (for (line <- fromFile(bamsIn).getLines) yield new File(line)).toList
    case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam and .list")
  }

  def getMapOfBAMsForSample(bams: List[File]): scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = {
    var m = scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]

    for (bam <- bams) {
      val bamSamples: List[String] = getSamplesInBAM(bam)
      for (s <- bamSamples) {
        if (!m.contains(s))
          m += s -> scala.collection.mutable.Set.empty[File]

        m(s) += bam
      }
    }

      return m
  }

  def findBAMsForSamples(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[File] = {
    var s = scala.collection.mutable.Set.empty[File]

    for (sample <- samples) {
      if (sampleToBams.contains(sample))
        s ++= sampleToBams(sample)
    }

    val l: List[File] = Nil
    return l ++ s
  }
}
