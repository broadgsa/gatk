package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource
import java.io.File
import net.sf.picard.util.IntervalList
import net.sf.samtools.SAMFileHeader
import collection.JavaConversions._
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocSortedSet, GenomeLocParser}
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, ScatterGatherableFunction, ScatterFunction}
import org.broadinstitute.sting.queue.function.{QFunction, InProcessFunction}
import org.broadinstitute.sting.queue.QException

/**
 * An interval scatter function.
 */
class IntervalScatterFunction extends ScatterFunction with InProcessFunction {
  var splitByContig = false

  private var referenceSequence: File = _
  private var intervalField: ArgumentSource = _
  private var intervals: List[String] = Nil

  def isScatterGatherable(originalFunction: ScatterGatherableFunction) = {
    if (originalFunction.isInstanceOf[CommandLineGATK]) {
      val gatk = originalFunction.asInstanceOf[CommandLineGATK]
      gatk.reference_sequence != null
    } else false
  }

  def setScatterGatherable(originalFunction: ScatterGatherableFunction) = {
    val gatk = originalFunction.asInstanceOf[CommandLineGATK]
    this.referenceSequence = gatk.reference_sequence
    this.intervals ++= gatk.intervalsString
    this.intervals ++= gatk.intervals.map(_.toString)
    this.intervalField = QFunction.findField(originalFunction.getClass, "intervals")
  }

  def initCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    cloneFunction.setFieldValue(this.intervalField, List(new File("scatter.intervals")))
  }

  def bindCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    val scatterPart = cloneFunction.getFieldValue(this.intervalField)
            .asInstanceOf[List[File]]
            .map(file => IOUtils.subDir(cloneFunction.commandDirectory, file))
    cloneFunction.setFieldValue(this.intervalField, scatterPart)
    this.scatterParts ++= scatterPart
  }

  def run() = {
    IntervalScatterFunction.scatter(this.intervals, this.scatterParts, this.referenceSequence, this.splitByContig)
  }
}

object IntervalScatterFunction {
  def scatter(intervals: List[String], scatterParts: List[File], reference: File, splitByContig: Boolean) = {
    val referenceSource = new ReferenceDataSource(reference)
    GenomeLocParser.setupRefContigOrdering(referenceSource.getReference)
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

    if (locs == null || locs.size == 0)
      throw new QException("Locs produced an empty interval list: " + intervals.mkString(", "))

    if (splitByContig) {
      var contig: String = null
      for (loc <- locs) {
        if (contig != loc.getContig && (fileIndex + 1) < scatterParts.size) {
          if (fileIndex >= 0)
            intervalList.write(scatterParts(fileIndex))
          fileIndex += 1
          contig = loc.getContig
          intervalList = new IntervalList(fileHeader)
        }
        locIndex += 1
        intervalList.add(toInterval(loc, locIndex))
      }
      intervalList.write(scatterParts(fileIndex))
      if ((fileIndex + 1) != scatterParts.size)
        throw new QException("Only able to write contigs into %d of %d files.".format(fileIndex + 1, scatterParts.size))
    } else {
      var locsPerFile = locs.size / scatterParts.size
      val locRemainder = locs.size % scatterParts.size

      // At the start, put an extra loc per file
      locsPerFile += 1
      var locsLeftFile = 0

      for (loc <- locs) {
        if (locsLeftFile == 0) {
          if (fileIndex >= 0)
            intervalList.write(scatterParts(fileIndex))

          fileIndex += 1
          intervalList = new IntervalList(fileHeader)

          // When we have put enough locs into each file,
          // reduce the number of locs per file back
          // to the original calculated value.
          if (fileIndex == locRemainder)
            locsPerFile -= 1
          locsLeftFile = locsPerFile
        }
        locsLeftFile -= 1
        locIndex += 1
        intervalList.add(toInterval(loc, locIndex))
      }
      intervalList.write(scatterParts(fileIndex))
      if ((fileIndex + 1) != scatterParts.size)
        throw new QException("Only able to write intervals into %d of %d files.".format(fileIndex + 1, scatterParts.size))
    }
  }

  private def toInterval(loc: GenomeLoc, locIndex: Int) =
    new net.sf.picard.util.Interval(loc.getContig, loc.getStart.toInt, loc.getStop.toInt, false, "interval_" + locIndex)
}
