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
  private var intervalsField: ArgumentSource = _
  private var intervalsStringField: ArgumentSource = _
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
    this.intervalsField = QFunction.findField(originalFunction.getClass, "intervals")
    this.intervalsStringField = QFunction.findField(originalFunction.getClass, "intervalsString")
  }

  def initCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    cloneFunction.setFieldValue(this.intervalsField, List(new File("scatter.intervals")))
    cloneFunction.setFieldValue(this.intervalsStringField, List.empty[String])
  }

  def bindCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    val scatterPart = cloneFunction.getFieldValue(this.intervalsField)
            .asInstanceOf[List[File]]
            .map(file => IOUtils.absolute(cloneFunction.commandDirectory, file))
    cloneFunction.setFieldValue(this.intervalsField, scatterPart)
    this.scatterParts ++= scatterPart
  }

  def run() = {
    IntervalUtils.scatterIntervalArguments(this.referenceSequence, this.intervals, this.scatterParts, this.splitByContig)
  }
}
