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

package org.broadinstitute.gatk.queue.extensions.gatk

import org.broadinstitute.gatk.utils.interval.IntervalUtils
import java.io.File
import org.broadinstitute.gatk.utils.io.IOUtils
import org.broadinstitute.gatk.queue.function.scattergather.{CloneFunction, ScatterFunction}
import org.broadinstitute.gatk.utils.commandline._

trait GATKScatterFunction extends ScatterFunction {
  /* The runtime field to set for specifying intervals. */
  private final val intervalsField = "intervals"
  private final val intervalsStringField = "intervalsString"
  private final val excludeIntervalsField = "excludeIntervals"
  private final val excludeIntervalsStringField = "excludeIntervalsString"
  private final val intervalsSetRuleField = "interval_set_rule"
  private final val intervalMergingField = "interval_merging"
  private final val intervalPaddingField = "interval_padding"

  @Output(doc="Scatter function outputs")
  var scatterOutputFiles: Seq[File] = Nil

  /** The original GATK function. */
  protected var originalGATK: CommandLineGATK = _

  /** Whether the last scatter job should also include any unmapped reads. */
  var includeUnmapped: Boolean = _

  override def init() {
    this.originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]
    // If intervals have been specified check if unmapped is included
    if (this.originalGATK.intervals.size + this.originalGATK.intervalsString.size > 0)
      this.includeUnmapped = this.originalGATK.intervalsString.exists(interval => IntervalUtils.isUnmapped(interval))
  }

  override def isScatterGatherable = {
    this.originalGATK.reference_sequence != null
  }

  override def initCloneInputs(cloneFunction: CloneFunction, index: Int) {
    cloneFunction.setFieldValue(this.intervalsField, Seq(new File("scatter.intervals")))
    if (index == this.scatterCount && this.includeUnmapped)
      cloneFunction.setFieldValue(this.intervalsStringField, Seq("unmapped"))
    else
      cloneFunction.setFieldValue(this.intervalsStringField, Seq.empty[String])

    cloneFunction.setFieldValue(this.intervalsSetRuleField, null)
    cloneFunction.setFieldValue(this.intervalMergingField, null)
    cloneFunction.setFieldValue(this.intervalPaddingField, None)
    cloneFunction.setFieldValue(this.excludeIntervalsField, Seq.empty[File])
    cloneFunction.setFieldValue(this.excludeIntervalsStringField, Seq.empty[String])
  }

  override def bindCloneInputs(cloneFunction: CloneFunction, index: Int) {
    val scatterPart = cloneFunction.getFieldValue(this.intervalsField)
            .asInstanceOf[Seq[File]]
            .map(file => IOUtils.absolute(cloneFunction.commandDirectory, file))
    cloneFunction.setFieldValue(this.intervalsField, scatterPart)
    this.scatterOutputFiles ++= scatterPart
  }

  /**
   * @return true if all interval files exist.
   */
  protected def intervalFilesExist = {
    !(this.originalGATK.intervals ++ this.originalGATK.excludeIntervals).exists(interval => !interval.exists())
  }

  /**
   * @return the maximum number of intervals or this.scatterCount if the maximum can't be determined ahead of time.
   */
  protected def maxIntervals: Int
}

object GATKScatterFunction {
  var gatkIntervalsCache = Seq.empty[GATKIntervals]

  def getGATKIntervals(originalFunction: CommandLineGATK) = {
    val gatkIntervals = new GATKIntervals(originalFunction)
    gatkIntervalsCache.find(_ == gatkIntervals) match {
      case Some(existingGatkIntervals) => existingGatkIntervals
      case None =>
        gatkIntervalsCache :+= gatkIntervals
        gatkIntervals
    }
  }
}
