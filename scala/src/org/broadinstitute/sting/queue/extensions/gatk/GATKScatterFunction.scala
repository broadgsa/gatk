/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.utils.interval.IntervalUtils
import java.io.File
import collection.JavaConversions._
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, ScatterGatherableFunction, ScatterFunction}
import org.broadinstitute.sting.commandline.{Output, ArgumentSource}

trait GATKScatterFunction extends ScatterFunction {
  /** The total number of clone jobs that will be created. */
  var scatterCount: Int = _

  /** The reference sequence for the GATK function. */
  protected var referenceSequence: File = _

  /** The runtime field to set for specifying an interval file. */
  protected var intervalsField: ArgumentSource = _

  /** The runtime field to set for specifying an interval string. */
  protected var intervalsStringField: ArgumentSource = _

  /** The list of interval files ("/path/to/interval.list") or interval strings ("chr1", "chr2") to parse into smaller parts. */
  protected var intervals: List[String] = Nil

  /** Whether the last scatter job should also include any unmapped reads. */
  protected var includeUnmapped: Boolean = _

  @Output(doc="Scatter function outputs")
  var scatterOutputFiles: List[File] = Nil

  /**
   * Checks if the function is scatter gatherable.
   * @param originalFunction Function to check.
   * @return true if the function is a GATK function with the reference sequence set.
   * @throws IllegalArgumentException if -BTI or -BTIMR are set.  QScripts should not try to scatter gather with those option set.
   */
  def isScatterGatherable(originalFunction: ScatterGatherableFunction): Boolean = {
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
    this.referenceSequence = gatk.reference_sequence
    if (gatk.intervals.isEmpty && gatk.intervalsString.isEmpty) {
      this.intervals ++= IntervalUtils.distinctContigs(this.referenceSequence).toList
    } else {
      this.intervals ++= gatk.intervals.map(_.toString)
      this.intervals ++= gatk.intervalsString.filterNot(interval => IntervalUtils.isUnmapped(interval))
      this.includeUnmapped = gatk.intervalsString.exists(interval => IntervalUtils.isUnmapped(interval))
    }

    this.scatterCount = originalFunction.scatterCount
    this.scatterCount = this.scatterCount min this.maxIntervals
  }

  def initCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    cloneFunction.setFieldValue(this.intervalsField, List(new File("scatter.intervals")))
    if (index == this.scatterCount && this.includeUnmapped)
      cloneFunction.setFieldValue(this.intervalsStringField, List("unmapped"))
    else
      cloneFunction.setFieldValue(this.intervalsStringField, List.empty[String])
  }

  def bindCloneInputs(cloneFunction: CloneFunction, index: Int) = {
    val scatterPart = cloneFunction.getFieldValue(this.intervalsField)
            .asInstanceOf[List[File]]
            .map(file => IOUtils.absolute(cloneFunction.commandDirectory, file))
    cloneFunction.setFieldValue(this.intervalsField, scatterPart)
    this.scatterOutputFiles ++= scatterPart
  }

  /**
   * Returns true if all interval files exist.
   */
  protected def intervalFilesExist = {
    !this.intervals.exists(interval => IntervalUtils.isIntervalFile(interval, false) && !new File(interval).exists)
  }

  /**
   * Returns the maximum number of intervals or this.scatterCount if the maximum can't be determined ahead of time.
   * @return the maximum number of intervals or this.scatterCount if the maximum can't be determined ahead of time.
   */
  protected def maxIntervals = this.scatterCount
}
