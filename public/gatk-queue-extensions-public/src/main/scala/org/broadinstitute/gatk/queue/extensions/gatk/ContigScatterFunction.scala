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

import collection.JavaConversions._
import org.broadinstitute.gatk.utils.interval.IntervalUtils
import org.broadinstitute.gatk.queue.function.InProcessFunction

/**
 * Splits intervals by contig instead of evenly.
 */
class ContigScatterFunction extends GATKScatterFunction with InProcessFunction {

  override def scatterCount = if (intervalFilesExist) super.scatterCount min this.maxIntervals else super.scatterCount

  protected override def maxIntervals = {
    GATKScatterFunction.getGATKIntervals(this.originalGATK).contigs.size
  }

  def run() {
    val gi = GATKScatterFunction.getGATKIntervals(this.originalGATK)
    IntervalUtils.scatterContigIntervals(gi.samFileHeader, gi.locs, this.scatterOutputFiles)
  }
}

