package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.function.scattergather.{ScatterGatherableFunction, ScatterFunction}
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource
import java.io.File
import net.sf.picard.util.IntervalList
import net.sf.samtools.SAMFileHeader
import collection.JavaConversions._
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocSortedSet, GenomeLocParser}

/**
 * An interval scatter function.
 */
class IntervalScatterFunction extends ScatterFunction with InProcessFunction {
  var splitByContig = false

  private var referenceSequence: File = _
  private var intervals: List[String] = Nil

  override def setOriginalFunction(originalFunction: ScatterGatherableFunction, scatterField: ArgumentSource) = {
    val gatk = originalFunction.asInstanceOf[CommandLineGATK]
    referenceSequence = gatk.reference_sequence
    intervals = gatk.intervalsString
    if (gatk.intervals != null)
      intervals ::= gatk.intervals.toString
  }

  def run() = {
    val referenceSource = new ReferenceDataSource(referenceSequence)
    GenomeLocParser.setupRefContigOrdering(referenceSource.getReference);
    val locs = {
      // TODO: Abstract genome analysis engine has richer logic for parsing.  We need to use it!
      if (intervals.size == 0) {
        GenomeLocSortedSet.createSetFromSequenceDictionary(referenceSource.getReference.getSequenceDictionary).toList
      } else {
        IntervalUtils.parseIntervalArguments(intervals, false)
      }
    }

    val fileHeader = new SAMFileHeader
    fileHeader.setSequenceDictionary(referenceSource.getReference.getSequenceDictionary)

    var intervalList: IntervalList = null
    var fileIndex = -1
    var locIndex = 0

    if (splitByContig) {
      var contig: String = null
      for (loc <- locs) {
        if (contig != loc.getContig && (fileIndex + 1) < scatterParts.size) {
          if (fileIndex >= 0)
            intervalList.write(scatterParts(fileIndex))
          fileIndex += 1
          intervalList = new IntervalList(fileHeader)
        }
        locIndex += 1
        intervalList.add(toInterval(loc, locIndex))
      }
      intervalList.write(scatterParts(fileIndex))
    } else {
      var locsPerFile = locs.size / this.scatterParts.size
      if (locs.size % this.scatterParts.size != 0) locsPerFile += 1
      for (loc <- locs) {
        if (locIndex % locsPerFile == 0) {
          if (fileIndex >= 0)
            intervalList.write(scatterParts(fileIndex))
          fileIndex += 1
          intervalList = new IntervalList(fileHeader)
        }
        locIndex += 1
        intervalList.add(toInterval(loc, locIndex))
      }
      intervalList.write(scatterParts(fileIndex))
    }
  }

  private def toInterval(loc: GenomeLoc, locIndex: Int) =
    new net.sf.picard.util.Interval(loc.getContig, loc.getStart.toInt, loc.getStop.toInt, true, "interval_" + locIndex)
}
