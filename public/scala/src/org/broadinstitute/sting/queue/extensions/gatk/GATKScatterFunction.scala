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
import org.broadinstitute.sting.utils.io.IOUtils
import org.broadinstitute.sting.queue.function.scattergather.{CloneFunction, ScatterFunction}
import org.broadinstitute.sting.commandline.Output

trait GATKScatterFunction extends ScatterFunction {
  /** The runtime field to set for specifying an interval file. */
  private final val intervalsField = "intervals"

  /** The runtime field to set for specifying an interval string. */
  private final val intervalsStringField = "intervalsString"

  @Output(doc="Scatter function outputs")
  var scatterOutputFiles: List[File] = Nil

  /** The original GATK function. */
  protected var originalGATK: CommandLineGATK = _

  /** The reference sequence for the GATK function. */
  protected var referenceSequence: File = _

  /** The list of interval files ("/path/to/interval.list") or interval strings ("chr1", "chr2") to parse into smaller parts. */
  protected var intervals: List[String] = Nil

  /** Whether the last scatter job should also include any unmapped reads. */
  protected var includeUnmapped: Boolean = _

  override def init() {
    this.originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]
    this.referenceSequence = this.originalGATK.reference_sequence
    if (this.originalGATK.intervals.isEmpty && (this.originalGATK.intervalsString == null || this.originalGATK.intervalsString.isEmpty)) {
      this.intervals ++= GATKScatterFunction.getGATKIntervals(this.referenceSequence, List.empty[String]).contigs
    } else {
      this.intervals ++= this.originalGATK.intervals.map(_.toString)
      this.intervals ++= this.originalGATK.intervalsString.filterNot(interval => IntervalUtils.isUnmapped(interval))
      this.includeUnmapped = this.originalGATK.intervalsString.exists(interval => IntervalUtils.isUnmapped(interval))
    }
  }

  override def isScatterGatherable = {
    this.originalGATK.reference_sequence != null
  }

  override def initCloneInputs(cloneFunction: CloneFunction, index: Int) {
    cloneFunction.setFieldValue(this.intervalsField, List(new File("scatter.intervals")))
    if (index == this.scatterCount && this.includeUnmapped)
      cloneFunction.setFieldValue(this.intervalsStringField, List("unmapped"))
    else
      cloneFunction.setFieldValue(this.intervalsStringField, List.empty[String])
  }

  override def bindCloneInputs(cloneFunction: CloneFunction, index: Int) {
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
  protected def maxIntervals: Int
}

object GATKScatterFunction {
  var gatkIntervals = List.empty[GATKIntervals]

  def getGATKIntervals(reference: File, intervals: List[String]) = {
    gatkIntervals.find(gi => gi.reference == reference && gi.intervals == intervals) match {
      case Some(gi) => gi
      case None =>
        val gi = new GATKIntervals(reference, intervals)
        gatkIntervals :+= gi
        gi
    }
  }
}
