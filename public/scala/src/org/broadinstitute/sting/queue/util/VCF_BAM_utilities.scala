package org.broadinstitute.sting.queue.util

import java.io.File
import org.apache.commons.io.FilenameUtils
import scala.io.Source._
import net.sf.samtools.SAMFileReader
import org.broadinstitute.sting.utils.codecs.vcf.{VCFHeader, VCFCodec}
import scala.collection.JavaConversions._
import org.broad.tribble.{FeatureCodec, AbstractFeatureReader}
import org.broadinstitute.sting.utils.variantcontext.VariantContext

object VCF_BAM_utilities {

  def getSamplesFromVCF(vcfFile: File): List[String] = {
    val codec: FeatureCodec[VariantContext] = new VCFCodec().asInstanceOf[FeatureCodec[VariantContext]]
    AbstractFeatureReader.getFeatureReader(vcfFile.getPath, codec).getHeader.asInstanceOf[VCFHeader].getGenotypeSamples.toList
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
