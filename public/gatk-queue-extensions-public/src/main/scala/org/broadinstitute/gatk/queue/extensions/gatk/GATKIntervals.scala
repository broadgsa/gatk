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

import java.io.File
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.interval.{IntervalSetRule, IntervalMergingRule, IntervalUtils}
import org.broadinstitute.gatk.engine.datasources.reference.ReferenceDataSource
import htsjdk.samtools.SAMFileHeader
import org.broadinstitute.gatk.utils.GenomeLoc
import org.broadinstitute.gatk.utils.commandline._
import htsjdk.tribble.Feature

case class GATKIntervals(reference: File, intervals: Seq[File], intervalsString: Seq[String],
                         intervalSetRule: IntervalSetRule, intervalMergingRule: IntervalMergingRule, intervalPadding: Option[Int],
                         excludeIntervals: Seq[File], excludeIntervalsString: Seq[String]) {

  def this(gatk: CommandLineGATK) = this(
    gatk.reference_sequence,
    gatk.intervals, gatk.intervalsString,
    gatk.interval_set_rule, gatk.interval_merging, gatk.interval_padding,
    gatk.excludeIntervals, gatk.excludeIntervalsString)

  private lazy val referenceDataSource = new ReferenceDataSource(reference)

  lazy val samFileHeader = {
    val header = new SAMFileHeader
    header.setSequenceDictionary(referenceDataSource.getReference.getSequenceDictionary)
    header
  }

  lazy val locs: java.util.List[GenomeLoc] = {
    val includeIntervalBindings = this.intervals.map(GATKIntervals.createBinding(_, "intervals")) ++
      this.intervalsString.map(GATKIntervals.createBinding(_, "intervalsString"))
    val excludeIntervalBindings = this.excludeIntervals.map(GATKIntervals.createBinding(_, "excludeIntervals")) ++
      this.excludeIntervalsString.map(GATKIntervals.createBinding(_, "excludeIntervalsString"))

    IntervalUtils.parseIntervalBindings(
      referenceDataSource.getReference,
      includeIntervalBindings,
      intervalSetRule, intervalMergingRule, intervalPadding.getOrElse(0),
      excludeIntervalBindings).toList
  }

  lazy val contigs = locs.map(_.getContig).distinct.toSeq
}

object GATKIntervals {
  def copyIntervalArguments(src: CommandLineGATK, dst: CommandLineGATK) {
    dst.reference_sequence = src.reference_sequence
    dst.intervals = src.intervals
    dst.intervalsString = src.intervalsString
    dst.interval_set_rule = src.interval_set_rule
    dst.interval_merging = src.interval_merging
    dst.interval_padding = src.interval_padding
    dst.excludeIntervals = src.excludeIntervals
    dst.excludeIntervalsString = src.excludeIntervalsString
  }

  private def createBinding(interval: File, argumentName: String): IntervalBinding[Feature] = {
    val tags = interval match {
      case taggedFile: TaggedFile => ParsingMethod.parseTags(argumentName, taggedFile.tag)
      case file: File => new Tags
    }
    createBinding(interval.getAbsolutePath, argumentName, tags)
  }

  private def createBinding(interval: String, argumentName: String): IntervalBinding[Feature] = {
    createBinding(interval, argumentName, new Tags)
  }

  private def createBinding(interval: String, argumentName: String, tags: Tags): IntervalBinding[Feature] = {
    ArgumentTypeDescriptor.parseBinding(new ArgumentMatchStringValue(interval), classOf[Feature], classOf[IntervalBinding[Feature]], argumentName, tags, argumentName).asInstanceOf[IntervalBinding[Feature]]
  }
}
