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
import org.testng.Assert
import org.testng.annotations.{DataProvider, Test}
import org.broadinstitute.gatk.utils.BaseTest
import org.broadinstitute.gatk.engine.datasources.reference.ReferenceDataSource
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile
import org.broadinstitute.gatk.utils.{GenomeLocSortedSet, GenomeLocParser}
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.interval.IntervalUtils
import org.broadinstitute.gatk.utils.exceptions.UserException

class GATKIntervalsUnitTest {
  private final lazy val hg18Reference = new File(BaseTest.hg18Reference)
  private final lazy val hg18GenomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(hg18Reference))
  private final lazy val hg18ReferenceLocs = GenomeLocSortedSet.
    createSetFromSequenceDictionary(new ReferenceDataSource(hg18Reference).getReference.getSequenceDictionary).toList
  private final lazy val hg19GenomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(hg19Reference))

  private final lazy val hg19Reference = new File(BaseTest.hg19Reference)

  @Test
  def testWithIntervals() {
    val chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1")
    val chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:2-3")
    val chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:3-5")

    val gi = createGATKIntervals(hg18Reference, Seq("chr1:1-1", "chr2:2-3", "chr3:3-5"))
    Assert.assertEquals(gi.locs.toSeq, Seq(chr1, chr2, chr3))
    Assert.assertEquals(gi.contigs, Seq("chr1", "chr2", "chr3"))
  }

  @Test(timeOut = 30000L)
  def testIntervalFile() {
    val gi = createGATKIntervals(hg19Reference, Seq(BaseTest.hg19Intervals))
    Assert.assertEquals(gi.locs.size, 189894)
    // Timeout check is because of bad:
    //   for(Item item: javaConvertedScalaList)
    // This for loop is actually an O(N^2) operation as the iterator calls the
    // O(N) javaConvertedScalaList.size() for each iteration of the loop.
    Assert.assertEquals(IntervalUtils.splitFixedIntervals(gi.locs, 189894).size(), 189894)
    Assert.assertEquals(gi.contigs.size, 24)
  }

  @Test
  def testEmptyIntervals() {
    val gi = createGATKIntervals(hg18Reference, Nil)
    Assert.assertEquals(gi.locs, hg18ReferenceLocs)
    Assert.assertEquals(gi.contigs.size, hg18ReferenceLocs.size)
  }

  @Test
  def testContigCounts() {
    Assert.assertEquals(createGATKIntervals(hg18Reference, Nil).contigs, hg18ReferenceLocs.map(_.getContig))
    Assert.assertEquals(createGATKIntervals(hg18Reference, Seq("chr1", "chr2", "chr3")).contigs, Seq("chr1", "chr2", "chr3"))
    Assert.assertEquals(createGATKIntervals(hg18Reference, Seq("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2")).contigs, Seq("chr1", "chr2", "chr3"))
  }

  @DataProvider(name="sortAndMergeIntervals")
  def getSortAndMergeIntervals: Array[Array[AnyRef]] = {
    Array(
      Array(Seq("chr1:1-10", "chr1:1-10", "chr1:1-10"), Seq("chr1:1-10")),
      Array(Seq("chr1:1-10", "chr1:1-11", "chr1:1-12"), Seq("chr1:1-12")),
      Array(Seq("chr1:1-10", "chr1:11-20", "chr1:21-30"), Seq("chr1:1-30")),
      Array(Seq("chr1:1-10", "chr1:10-20", "chr1:21-30"), Seq("chr1:1-30")),
      Array(Seq("chr1:1-9", "chr1:21-30", "chr1:11-20"), Seq("chr1:1-9", "chr1:11-30"))
    ).asInstanceOf[Array[Array[AnyRef]]]
  }

  @Test(dataProvider="sortAndMergeIntervals")
  def testSortAndMergeIntervals(unmerged: Seq[String], expected: Seq[String]) {
    Assert.assertEquals(createGATKIntervals(hg18Reference, unmerged).locs.toSeq, expected.map(hg18GenomeLocParser.parseGenomeLoc(_)))
  }

  @DataProvider(name="taggedFiles")
  def getTaggedFiles: Array[Array[AnyRef]] = {
    Array(
      Array(hg18Reference, BaseTest.privateTestDir + "small_unmerged_gatk_intervals.list", null, Seq("chr1:1-10")),
      Array(hg18Reference, BaseTest.privateTestDir + "small_unmerged_gatk_intervals.list", "", Seq("chr1:1-10")),
      Array(hg18Reference, BaseTest.privateTestDir + "small_unmerged_gatk_intervals.list", "myList", Seq("chr1:1-10")),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", null, Seq("1:897475-897481", "1:10001292")),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "", Seq("1:897475-897481", "1:10001292")),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "myVcf", Seq("1:897475-897481", "1:10001292")),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "VCF", Seq("1:897475-897481", "1:10001292")),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "myVcf,VCF", Seq("1:897475-897481", "1:10001292")),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", null, Seq("20:1-999", "20:1002-2000", "22:1001-6000")),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "", Seq("20:1-999", "20:1002-2000", "22:1001-6000")),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "myBed", Seq("20:1-999", "20:1002-2000", "22:1001-6000")),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "BED", Seq("20:1-999", "20:1002-2000", "22:1001-6000")),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "myBed,BED", Seq("20:1-999", "20:1002-2000", "22:1001-6000"))
    )
  }

  @Test(dataProvider="taggedFiles")
  def testTaggedFiles(reference: File, file: String, tags: String, expected: Seq[String]) {
    val gatk = new CommandLineGATK
    gatk.reference_sequence = reference
    gatk.intervals = Seq(new TaggedFile(file, tags))
    val parser = if (reference == hg18Reference) hg18GenomeLocParser else hg19GenomeLocParser
    Assert.assertEquals(new GATKIntervals(gatk).locs.toSeq, expected.map(parser.parseGenomeLoc(_)))
  }

  @DataProvider(name="badTaggedFiles")
  def getBadTaggedFiles: Array[Array[AnyRef]] = {
    Array(
      Array(hg18Reference, BaseTest.privateTestDir + "small_unmerged_gatk_intervals.list", "VCF"),
      Array(hg18Reference, BaseTest.privateTestDir + "small_unmerged_gatk_intervals.list", "too,many,tags"),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "BED"),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "VCF,myVCF"),
      Array(hg19Reference, BaseTest.privateTestDir + "small.indel.test.vcf", "myVCF,VCF,extra"),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "VCF"),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "BED,myBed"),
      Array(hg19Reference, BaseTest.privateTestDir + "sampleBedFile.bed", "myBed,BED,extra")
    ).asInstanceOf[Array[Array[AnyRef]]]
  }

  @Test(dataProvider = "badTaggedFiles", expectedExceptions = Array(classOf[UserException]))
  def testBadTaggedFiles(reference: File, file: String, tags: String) {
    testTaggedFiles(reference, file, tags, Nil)
  }

  private def createGATKIntervals(reference: File, intervals: Seq[String]) = {
    val gatk = new CommandLineGATK
    gatk.reference_sequence = reference
    gatk.intervalsString = intervals
    new GATKIntervals(gatk)
  }
}
