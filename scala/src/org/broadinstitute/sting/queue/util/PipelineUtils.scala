package org.broadinstitute.sting.queue.util

import net.sf.picard.reference.ReferenceSequenceFileFactory
import java.io.File
import org.broadinstitute.sting.utils.GenomeLocParser
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils

class PipelineUtils {

}

object PipelineUtils{

  def smartSplitContigs(reference: File, intervals: File, sets: Int) : List[List[String]] = {
    GenomeLocParser.setupRefContigOrdering(ReferenceSequenceFileFactory.getReferenceSequenceFile(reference))
    val targets = IntervalUtils.parseIntervalArguments(List(intervals.getAbsolutePath), false)

    // Build up a map of contigs with sizes.
    var contigSizes = Map.empty[String, Long]
    // todo -- make this look like functional code for Q's sake
    //targets.foreach( loc => { contigSizes += loc -> { contigSizes.get(loc.getContig) match { case Some(size) => size + loc.size case None => loc.size } } })

    for (loc <- targets) {
      val contig = loc.getContig
      val contigSize = loc.size
      contigSizes += contig -> {
        contigSizes.get(contig) match {
          case Some(size) => size + contigSize
          case None => contigSize
        }
      }
    }

    // Keep a list of pairs of sizes with lists of contigs.
    var splitContigs = List.empty[(Long, List[String])]
    for ((contigName, contigSize) <- contigSizes) {
      if (splitContigs.size < sets) {
        // If there are fewer than the requested number of sets, just add this contig.
        splitContigs :+= contigSize -> List(contigName)
      } else {
        // If there is already a number of sets
        // sort the contigs to get the smallest one first.
        splitContigs = splitContigs.sortBy{case (size, contigs) => size}
        // Update the pair with the new contig size and name.
        var smallContigs = splitContigs.head
        smallContigs = (smallContigs._1 + contigSize) -> (smallContigs._2 :+ contigName)
        // Re add the pair to the list.
        splitContigs = smallContigs :: splitContigs.tail
      }
    }

    splitContigs.map{case (size, contigs) => contigs}
  }
}