package org.broadinstitute.sting.queue.util

import tools.nsc.io.Process
import java.io.File
import org.apache.commons.io.FilenameUtils
import scala.io.Source._

object BAMutilities {

  def getSamplesFromVCF(calls: File): List[String] = {
    val getSamplesCommand: String = "cat " + calls.getPath() + " | grep '^#CHROM' | head -1 | awk '{for (i = 10; i <= NF; i++) print $i}'"
    //println("getSamplesCommand: " + getSamplesCommand)
    return Process(getSamplesCommand).iterator toList
  }

  def parseBamsInput(bamsIn: File): List[File] = FilenameUtils.getExtension(bamsIn.getPath) match {
    case "bam" => return List(bamsIn)
    case "list" => return (for (line <- fromFile(bamsIn).getLines) yield new File(line)).toList
    case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam and .list")
  }

  def getMapOfBamsForSample(bams: List[File]): scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = bams match {
    case Nil => return scala.collection.mutable.Map.empty[String, scala.collection.mutable.Set[File]]

    case x :: y =>
      val m: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = getMapOfBamsForSample(y)

      val getBamSamplesCommand: String = "samtools view -H " + x.getPath() + " | grep '^@RG' | awk '{for (i = 1; i <= NF; i++) if (substr($i,1,3) == \"SM:\") print substr($i,4)}' | sort | uniq"
      //println("getBamSamplesCommand: " + getBamSamplesCommand)
      val bamSamples: List[String] = Process(getBamSamplesCommand).iterator toList

      for (s <- bamSamples) {
        if (!m.contains(s))
          m += s -> scala.collection.mutable.Set.empty[File]

        m(s) = m(s) + x
      }

      return m
  }

  def findBamsForSamples(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[File] = {
    val l: List[File] = Nil
    l ++ findBamsForSamplesHelper(samples, sampleToBams)
  }

  def findBamsForSamplesHelper(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): scala.collection.mutable.Set[File] = samples match {
    case Nil => scala.collection.mutable.Set.empty[File]

    case x :: y =>
      var bamsForSampleX: scala.collection.mutable.Set[File] = scala.collection.mutable.Set.empty[File]
      if (sampleToBams.contains(x))
        bamsForSampleX = sampleToBams(x)
      return bamsForSampleX ++ findBamsForSamplesHelper(y, sampleToBams)
  }
}