package org.broadinstitute.sting.queue.extensions.gatk

import collection.JavaConversions._
import java.io.File
import org.testng.Assert
import org.broadinstitute.sting.BaseTest
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.queue.QException
import net.sf.picard.reference.IndexedFastaSequenceFile
import org.testng.annotations.{Test, BeforeMethod}
import org.broadinstitute.sting.utils.GenomeLocParser

class IntervalScatterFunctionUnitTest extends BaseTest {
  private def reference = new File(BaseTest.b36KGReference)
  private var header: IndexedFastaSequenceFile = _
  private var genomeLocParser: GenomeLocParser = _

  @BeforeMethod
  def setup() {
    header = new IndexedFastaSequenceFile(reference)
    genomeLocParser = new GenomeLocParser(header.getSequenceDictionary())
  }

  @Test
  def testCountContigs = {
    Assert.assertEquals(List("1"), IntervalScatterFunction.distinctContigs(reference, List(BaseTest.validationDataLocation + "chr1_b36_pilot3.interval_list")))
    Assert.assertEquals(List("1","2","3"), IntervalScatterFunction.distinctContigs(reference, List("1:1-1", "2:1-1", "3:2-2")))
    Assert.assertEquals(List("1","2","3"), IntervalScatterFunction.distinctContigs(reference, List("2:1-1", "1:1-1", "3:2-2")))
  }

  @Test
  def testBasicScatter = {
    val chr1 = genomeLocParser.parseGenomeInterval("1")
    val chr2 = genomeLocParser.parseGenomeInterval("2")
    val chr3 = genomeLocParser.parseGenomeInterval("3")

    val files = (1 to 3).toList.map(index => new File(testDir + "basic." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, List("1", "2", "3"), files, false)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testScatterLessFiles = {
    val chr1 = genomeLocParser.parseGenomeInterval("1")
    val chr2 = genomeLocParser.parseGenomeInterval("2")
    val chr3 = genomeLocParser.parseGenomeInterval("3")
    val chr4 = genomeLocParser.parseGenomeInterval("4")

    val files = (1 to 3).toList.map(index => new File(testDir + "less." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, List("1", "2", "3", "4"), files, false)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(2, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs1.get(1))
    Assert.assertEquals(chr3, locs2.get(0))
    Assert.assertEquals(chr4, locs3.get(0))
  }

  @Test(expectedExceptions=Array(classOf[QException]))
  def testScatterMoreFiles = {
    val files = (1 to 3).toList.map(index => new File(testDir + "more." + index + ".intervals"))
    IntervalScatterFunction.scatter(reference, List("1", "2"), files, false)
  }

  @Test
  def testScatterIntervals = {
    val intervals = List("1:1-2", "1:4-5", "2:1-1", "3:2-2")
    val chr1a = genomeLocParser.parseGenomeInterval("1:1-2")
    val chr1b = genomeLocParser.parseGenomeInterval("1:4-5")
    val chr2 = genomeLocParser.parseGenomeInterval("2:1-1")
    val chr3 = genomeLocParser.parseGenomeInterval("3:2-2")

    val files = (1 to 3).toList.map(index => new File(testDir + "split." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, intervals, files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(2, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1a, locs1.get(0))
    Assert.assertEquals(chr1b, locs1.get(1))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testScatterOrder = {
    val intervals = List("2:1-1", "1:1-1", "3:2-2")
    val chr1 = genomeLocParser.parseGenomeInterval("1:1-1")
    val chr2 = genomeLocParser.parseGenomeInterval("2:1-1")
    val chr3 = genomeLocParser.parseGenomeInterval("3:2-2")

    val files = (1 to 3).toList.map(index => new File(testDir + "split." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, intervals, files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testBasicScatterByContig = {
    val chr1 = genomeLocParser.parseGenomeInterval("1")
    val chr2 = genomeLocParser.parseGenomeInterval("2")
    val chr3 = genomeLocParser.parseGenomeInterval("3")

    val files = (1 to 3).toList.map(index => new File(testDir + "contig_basic." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, List("1", "2", "3"), files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testScatterByContigLessFiles = {
    val chr1 = genomeLocParser.parseGenomeInterval("1")
    val chr2 = genomeLocParser.parseGenomeInterval("2")
    val chr3 = genomeLocParser.parseGenomeInterval("3")
    val chr4 = genomeLocParser.parseGenomeInterval("4")

    val files = (1 to 3).toList.map(index => new File(testDir + "contig_less." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, List("1", "2", "3", "4"), files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(2, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
    Assert.assertEquals(chr4, locs3.get(1))
  }

  @Test(expectedExceptions=Array(classOf[QException]))
  def testScatterByContigMoreFiles = {
    val files = (1 to 3).toList.map(index => new File(testDir + "contig_more." + index + ".intervals"))
    IntervalScatterFunction.scatter(reference, List("1", "2"), files, true)
  }

  @Test
  def testScatterByContigIntervalsStart = {
    val intervals = List("1:1-2", "1:4-5", "2:1-1", "3:2-2")
    val chr1a = genomeLocParser.parseGenomeInterval("1:1-2")
    val chr1b = genomeLocParser.parseGenomeInterval("1:4-5")
    val chr2 = genomeLocParser.parseGenomeInterval("2:1-1")
    val chr3 = genomeLocParser.parseGenomeInterval("3:2-2")

    val files = (1 to 3).toList.map(index => new File(testDir + "contig_split_start." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, intervals, files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(2, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1a, locs1.get(0))
    Assert.assertEquals(chr1b, locs1.get(1))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testScatterByContigIntervalsMiddle = {
    val intervals = List("1:1-1", "2:1-2", "2:4-5", "3:2-2")
    val chr1 = genomeLocParser.parseGenomeInterval("1:1-1")
    val chr2a = genomeLocParser.parseGenomeInterval("2:1-2")
    val chr2b = genomeLocParser.parseGenomeInterval("2:4-5")
    val chr3 = genomeLocParser.parseGenomeInterval("3:2-2")

    val files = (1 to 3).toList.map(index => new File(testDir + "contig_split_middle." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, intervals, files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(2, locs2.size)
    Assert.assertEquals(1, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2a, locs2.get(0))
    Assert.assertEquals(chr2b, locs2.get(1))
    Assert.assertEquals(chr3, locs3.get(0))
  }

  @Test
  def testScatterByContigIntervalsEnd = {
    val intervals = List("1:1-1", "2:2-2", "3:1-2", "3:4-5")
    val chr1 = genomeLocParser.parseGenomeInterval("1:1-1")
    val chr2 = genomeLocParser.parseGenomeInterval("2:2-2")
    val chr3a = genomeLocParser.parseGenomeInterval("3:1-2")
    val chr3b = genomeLocParser.parseGenomeInterval("3:4-5")

    val files = (1 to 3).toList.map(index => new File(testDir + "contig_split_end." + index + ".intervals"))

    IntervalScatterFunction.scatter(reference, intervals, files, true)

    val locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(0).toString), false)
    val locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(1).toString), false)
    val locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser,List(files(2).toString), false)

    Assert.assertEquals(1, locs1.size)
    Assert.assertEquals(1, locs2.size)
    Assert.assertEquals(2, locs3.size)

    Assert.assertEquals(chr1, locs1.get(0))
    Assert.assertEquals(chr2, locs2.get(0))
    Assert.assertEquals(chr3a, locs3.get(0))
    Assert.assertEquals(chr3b, locs3.get(1))
  }
}
