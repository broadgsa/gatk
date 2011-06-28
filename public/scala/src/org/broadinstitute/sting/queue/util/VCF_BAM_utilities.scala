package org.broadinstitute.sting.queue.util

import java.io.File
import org.apache.commons.io.FilenameUtils
import scala.io.Source._
import net.sf.samtools.SAMFileReader
import org.broad.tribble.source.BasicFeatureSource
import org.broadinstitute.sting.utils.codecs.vcf.{VCFHeader, VCFCodec}
import scala.collection.JavaConversions._

object VCF_BAM_utilities {

  def getSamplesFromVCF(vcfFile: File): List[String] = {
    return BasicFeatureSource.getFeatureSource(vcfFile.getPath(), new VCFCodec()).getHeader().asInstanceOf[VCFHeader].getGenotypeSamples().toList
  }

  def getSamplesInBAM(bam: File): List[String] = {
    return new SAMFileReader(bam).getFileHeader().getReadGroups().toList.map(srgr => srgr.getSample()).toSet.toList
  }

  def parseBAMsInput(bamsIn: File): List[File] = FilenameUtils.getExtension(bamsIn.getPath) match {
    case "bam" => return List(bamsIn)
    case "list" => return (for (line <- fromFile(bamsIn).getLines) yield new File(line)).toList
    case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam and .list")
  }

  def getMapOfBAMsForSample(bams: List[File]): scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = bams match {
    case Nil => return scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]

    case x :: y =>
      val m: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = getMapOfBAMsForSample(y)
      val bamSamples: List[String] = getSamplesInBAM(x)

      for (s <- bamSamples) {
        if (!m.contains(s))
          m += s -> scala.collection.mutable.Set.empty[File]

        m(s) = m(s) + x
      }

      return m
  }

  def findBAMsForSamples(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[File] = {

    def findBAMsForSamplesHelper(samples: List[String]): scala.collection.mutable.Set[File] = samples match {
      case Nil => scala.collection.mutable.Set.empty[File]

      case x :: y =>
        var bamsForSampleX: scala.collection.mutable.Set[File] = scala.collection.mutable.Set.empty[File]
        if (sampleToBams.contains(x))
          bamsForSampleX = sampleToBams(x)
        return bamsForSampleX ++ findBAMsForSamplesHelper(y)
    }

    val l: List[File] = Nil
    return l ++ findBAMsForSamplesHelper(samples)
  }
}
