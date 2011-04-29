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
import org.testng.Assert
import org.testng.annotations.Test
import org.broadinstitute.sting.BaseTest
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile
import org.broadinstitute.sting.utils.{GenomeLocSortedSet, GenomeLocParser}
import collection.JavaConversions._

class GATKIntervalsUnitTest {
  private final lazy val reference = new File(BaseTest.hg18Reference)
  private final lazy val genomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(reference))
  private final lazy val referenceLocs = GenomeLocSortedSet.
    createSetFromSequenceDictionary(new ReferenceDataSource(reference).getReference.getSequenceDictionary).toList.toList

  @Test
  def testWithIntervals() {
    val chr1 = genomeLocParser.parseGenomeInterval("chr1:1-1")
    val chr2 = genomeLocParser.parseGenomeInterval("chr2:2-3")
    val chr3 = genomeLocParser.parseGenomeInterval("chr3:3-5")

    val gi = new GATKIntervals(reference, List("chr1:1-1", "chr2:2-3", "chr3:3-5"))
    Assert.assertEquals(gi.locs, List(chr1, chr2, chr3))
    Assert.assertEquals(gi.contigs, List("chr1", "chr2", "chr3"))
    Assert.assertEquals(gi.getSplits(2).toList, List(2, 3))
    Assert.assertEquals(gi.getSplits(3).toList, List(1, 2, 3))
  }

  @Test
  def testEmptyIntervals() {
    val gi = new GATKIntervals(reference, Nil)
    Assert.assertEquals(gi.locs, referenceLocs)
    Assert.assertEquals(gi.contigs.size, referenceLocs.size)
    Assert.assertEquals(gi.getSplits(2).toList, List(10, 45))
    Assert.assertEquals(gi.getSplits(4).toList, List(5, 10, 16, 45))
  }
}
