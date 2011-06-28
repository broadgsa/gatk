package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.BaseTest
import java.io.File
import org.broadinstitute.sting.utils.exceptions.UserException
import org.testng.Assert
import org.testng.annotations.Test

class IOUtilsUnitTest extends BaseTest {
  @Test
  def testGoodTempDir = {
    IOUtils.checkTempDir(new File("/tmp/queue"))
  }

  @Test(expectedExceptions=Array(classOf[UserException.BadTmpDir]))
  def testBadTempDir = {
    IOUtils.checkTempDir(new File("/tmp"))
  }

  @Test
  def testAbsoluteSubDir = {
    var subDir = IOUtils.absolute(new File("."), new File("/path/to/file"))
    Assert.assertEquals(subDir, new File("/path/to/file"))

    subDir = IOUtils.absolute(new File("/different/path"), new File("/path/to/file"))
    Assert.assertEquals(subDir, new File("/path/to/file"))

    subDir = IOUtils.absolute(new File("/different/path"), new File("."))
    Assert.assertEquals(subDir, new File("/different/path"))
  }

  @Test
  def testRelativeSubDir = {
    var subDir = IOUtils.absolute(new File("."), new File("path/to/file"))
    Assert.assertEquals(subDir.getCanonicalFile, new File("path/to/file").getCanonicalFile)

    subDir = IOUtils.absolute(new File("/different/path"), new File("path/to/file"))
    Assert.assertEquals(subDir, new File("/different/path/path/to/file"))
  }

  @Test
  def testDottedSubDir = {
    var subDir = IOUtils.absolute(new File("."), new File("path/../to/file"))
    Assert.assertEquals(subDir.getCanonicalFile, new File("path/../to/./file").getCanonicalFile)

    subDir = IOUtils.absolute(new File("."), new File("/path/../to/file"))
    Assert.assertEquals(subDir, new File("/path/../to/file"))

    subDir = IOUtils.absolute(new File("/different/../path"), new File("path/to/file"))
    Assert.assertEquals(subDir, new File("/different/../path/path/to/file"))

    subDir = IOUtils.absolute(new File("/different/./path"), new File("/path/../to/file"))
    Assert.assertEquals(subDir, new File("/path/../to/file"))
  }

  @Test
  def testTempDir = {
    val tempDir = IOUtils.tempDir("Q-Unit-Test", "", new File("queueTempDirToDelete"))
    Assert.assertTrue(tempDir.exists)
    Assert.assertFalse(tempDir.isFile)
    Assert.assertTrue(tempDir.isDirectory)
    val deleted = IOUtils.tryDelete(tempDir)
    Assert.assertTrue(deleted)
    Assert.assertFalse(tempDir.exists)
  }

  @Test
  def testDirLevel = {
    var dir = IOUtils.dirLevel(new File("/path/to/directory"), 1)
    Assert.assertEquals(dir, new File("/path"))

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 2)
    Assert.assertEquals(dir, new File("/path/to"))

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 3)
    Assert.assertEquals(dir, new File("/path/to/directory"))

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 4)
    Assert.assertEquals(dir, new File("/path/to/directory"))
  }

  @Test
  def testAbsolute = {
    var dir = IOUtils.absolute(new File("/path/./to/./directory/."))
    Assert.assertEquals(dir, new File("/path/to/directory"))

    dir = IOUtils.absolute(new File("/"))
    Assert.assertEquals(dir, new File("/"))

    dir = IOUtils.absolute(new File("/."))
    Assert.assertEquals(dir, new File("/"))

    dir = IOUtils.absolute(new File("/././."))
    Assert.assertEquals(dir, new File("/"))

    dir = IOUtils.absolute(new File("/./directory/."))
    Assert.assertEquals(dir, new File("/directory"))

    dir = IOUtils.absolute(new File("/./directory/./"))
    Assert.assertEquals(dir, new File("/directory"))

    dir = IOUtils.absolute(new File("/./directory./"))
    Assert.assertEquals(dir, new File("/directory."))

    dir = IOUtils.absolute(new File("/./.directory/"))
    Assert.assertEquals(dir, new File("/.directory"))
  }

  @Test
  def testTail = {
    val lines = List(
      "chr18_random	4262	3154410390	50	51",
      "chr19_random	301858	3154414752	50	51",
      "chr21_random	1679693	3154722662	50	51",
      "chr22_random	257318	3156435963	50	51",
      "chrX_random	1719168	3156698441	50	51")
    val tail = IOUtils.tail(new File(BaseTest.hg18Reference + ".fai"), 5)
    Assert.assertEquals(tail.size, 5)
    for (i <- 0 until 5)
      Assert.assertEquals(tail(i), lines(i))
  }
}
