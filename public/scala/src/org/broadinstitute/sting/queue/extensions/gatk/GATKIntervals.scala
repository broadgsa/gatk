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

import java.io.File
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource
import net.sf.samtools.SAMFileHeader
import java.util.Collections
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocSortedSet, GenomeLocParser}

case class GATKIntervals(reference: File, intervals: List[String]) {
  private lazy val referenceDataSource = new ReferenceDataSource(reference)
//  private var splitsBySize = Map.empty[Int, java.util.List[java.lang.Integer]]

  lazy val samFileHeader = {
    val header = new SAMFileHeader
    header.setSequenceDictionary(referenceDataSource.getReference.getSequenceDictionary)
    header
  }

  lazy val locs: java.util.List[GenomeLoc] = {
    val parser = new GenomeLocParser(referenceDataSource.getReference)
    val parsedLocs =
      if (intervals.isEmpty)
        GenomeLocSortedSet.createSetFromSequenceDictionary(samFileHeader.getSequenceDictionary).toList
      else
        IntervalUtils.parseIntervalArguments(parser, intervals)
    Collections.sort(parsedLocs)
    Collections.unmodifiableList(parsedLocs)
  }

  lazy val contigs = locs.map(_.getContig).distinct.toList

//  def getSplits(size: Int) = {
//    splitsBySize.getOrElse(size, {
//      val splits: java.util.List[java.lang.Integer] = IntervalUtils.splitFixedIntervals(locs, size)
//      splitsBySize += size -> splits
//      splits
//    })
//  }
}
