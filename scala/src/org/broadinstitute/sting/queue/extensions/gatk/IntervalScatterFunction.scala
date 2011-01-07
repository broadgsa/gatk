package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.utils.interval.IntervalUtils
import java.io.File
import collection.JavaConversions._
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, ScatterGatherableFunction, ScatterFunction}
import org.broadinstitute.sting.queue.function.{QFunction, InProcessFunction}

/**
 * An interval scatter function.
 */
class IntervalScatterFunction extends ScatterFunction with InProcessFunction {
  var splitByContig = false

  /** The total number of clone jobs that will be created. */
  private var scatterCount: Int = _

  /** The reference sequence for the GATK function. */
  private var referenceSequence: File = _

  /** The runtime field to set for specifying an interval file. */
  private var intervalsField: ArgumentSource = _

  /** The runtime field to set for specifying an interval string. */
  private var intervalsStringField: ArgumentSource = _

  /** The list of interval files ("/path/to/interval.list") or interval strings ("chr1", "chr2") to parse into smaller parts. */
  private var intervals: List[String] = Nil

  /** Whether the laster scatter job should also include any unmapped reads. */
  private var includeUnmapped: Boolean = _

  /**
   * Checks if the function is scatter gatherable.
   * @param originalFunction Function to check.
   * @return true if the function is a GATK function with the reference sequence set.
   * @throws IllegalArgumentException if -BTI or -BTIMR are set.  QScripts should not try to scatter gather with those option set.
   */
  def isScatterGatherable(originalFunction: ScatterGatherableFunction) = {
    if (originalFunction.isInstanceOf[CommandLineGATK]) {
      val gatk = originalFunction.asInstanceOf[CommandLineGATK]
      if ( gatk.BTI != null && gatk.BTIMR == null) throw new IllegalArgumentException("BTI requires BTIMR for use with scatter-gather (recommended: INTERSECTION)")
      gatk.reference_sequence != null
    } else false
  }

  /**
   * Sets the scatter gatherable function.
   * @param originalFunction Function to bind.
   */
  def setScatterGatherable(originalFunction: ScatterGatherableFunction) = {
    val gatk = originalFunction.asInstanceOf[CommandLineGATK]
    this.intervalsField = QFunction.findField(originalFunction.getClass, "intervals")
    this.intervalsStringField = QFunction.findField(originalFunction.getClass, "intervalsString")
    this.scatterCount = originalFunction.scatterCount
    this.referenceSequence = gatk.reference_sequence
    if (gatk.intervals.isEmpty && gatk.intervalsString.isEmpty) {
      this.intervals ++= IntervalUtils.distinctContigs(this.referenceSequence).toList
      this.includeUnmapped = this.splitByContig
    } else {
      this.intervals ++= gatk.intervals.map(_.toString)
      this.intervals ++= gatk.intervalsString.filterNot(interval => IntervalUtils.isUnmapped(interval))
      this.includeUnmapped = gatk.intervalsString.exists(interval => IntervalUtils.isUnmapped(interval))
    }
  }

  def initCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    cloneFunction.setFieldValue(this.intervalsField, List(new File("scatter.intervals")))
    if (index == scatterCount && includeUnmapped)
      cloneFunction.setFieldValue(this.intervalsStringField, List("unmapped"))
    else
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
