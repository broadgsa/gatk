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
  private final lazy val hg18Reference = new File(BaseTest.hg18Reference)
  private final lazy val hg18GenomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(hg18Reference))
  private final lazy val hg18ReferenceLocs = GenomeLocSortedSet.
    createSetFromSequenceDictionary(new ReferenceDataSource(hg18Reference).getReference.getSequenceDictionary).toList

  private final lazy val hg19Reference = new File(BaseTest.hg19Reference)
  private final lazy val hg19GenomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(hg19Reference))
  private final lazy val hg19ReferenceLocs = GenomeLocSortedSet.
    createSetFromSequenceDictionary(new ReferenceDataSource(hg19Reference).getReference.getSequenceDictionary).toList

  @Test
  def testWithIntervals() {
    val chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1")
    val chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:2-3")
    val chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:3-5")

    val gi = new GATKIntervals(hg18Reference, List("chr1:1-1", "chr2:2-3", "chr3:3-5"))
    Assert.assertEquals(gi.locs.toList, List(chr1, chr2, chr3))
    Assert.assertEquals(gi.contigs, List("chr1", "chr2", "chr3"))
//    Assert.assertEquals(gi.getSplits(2).toList, List(2, 3))
//    Assert.assertEquals(gi.getSplits(3).toList, List(1, 2, 3))
  }

  @Test(timeOut = 30000)
  def testIntervalFile() {
    var gi = new GATKIntervals(hg19Reference, List(BaseTest.hg19Intervals))
    Assert.assertEquals(gi.locs.size, 189894)
    // Timeout check is because of bad:
    //   for(Item item: javaConvertedScalaList)
    // This for loop is actually an O(N^2) operation as the iterator calls the
    // O(N) javaConvertedScalaList.size() for each iteration of the loop.
    //Assert.assertEquals(gi.getSplits(gi.locs.size).size, 189894)
    Assert.assertEquals(gi.contigs.size, 24)
  }

  @Test
  def testEmptyIntervals() {
    val gi = new GATKIntervals(hg18Reference, Nil)
    Assert.assertEquals(gi.locs, hg18ReferenceLocs)
    Assert.assertEquals(gi.contigs.size, hg18ReferenceLocs.size)
//    Assert.assertEquals(gi.getSplits(2).toList, List(10, 45))
//    Assert.assertEquals(gi.getSplits(4).toList, List(5, 10, 16, 45))
  }

  @Test
  def testContigCounts() {
    Assert.assertEquals(new GATKIntervals(hg18Reference, Nil).contigs, hg18ReferenceLocs.map(_.getContig))
    Assert.assertEquals(new GATKIntervals(hg18Reference, List("chr1", "chr2", "chr3")).contigs, List("chr1", "chr2", "chr3"))
    Assert.assertEquals(new GATKIntervals(hg18Reference, List("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2")).contigs, List("chr1", "chr2", "chr3"))
  }
}
