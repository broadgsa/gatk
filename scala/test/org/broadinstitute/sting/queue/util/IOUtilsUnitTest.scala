package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.BaseTest
import org.junit.{Assert, Test}
import java.io.File

class IOUtilsUnitTest extends BaseTest {
  @Test
  def testAbsoluteSubDir = {
    var subDir = IOUtils.subDir(IOUtils.CURRENT_DIR, new File("/path/to/file"))
    Assert.assertEquals(new File("/path/to/file"), subDir)

    subDir = IOUtils.subDir(new File("/different/path"), new File("/path/to/file"))
    Assert.assertEquals(new File("/path/to/file"), subDir)

    subDir = IOUtils.subDir(new File("/different/path"), IOUtils.CURRENT_DIR)
    Assert.assertEquals(new File("/different/path"), subDir)
  }

  @Test
  def testRelativeSubDir = {
    var subDir = IOUtils.subDir(IOUtils.CURRENT_DIR, new File("path/to/file"))
    Assert.assertEquals(new File("path/to/file").getCanonicalFile, subDir.getCanonicalFile)

    subDir = IOUtils.subDir(new File("/different/path"), new File("path/to/file"))
    Assert.assertEquals(new File("/different/path/path/to/file"), subDir)
  }

  @Test
  def testDottedSubDir = {
    var subDir = IOUtils.subDir(IOUtils.CURRENT_DIR, new File("path/../to/file"))
    Assert.assertEquals(new File("path/../to/./file").getCanonicalFile, subDir.getCanonicalFile)

    subDir = IOUtils.subDir(IOUtils.CURRENT_DIR, new File("/path/../to/file"))
    Assert.assertEquals(new File("/path/../to/file"), subDir)

    subDir = IOUtils.subDir(new File("/different/../path"), new File("path/to/file"))
    Assert.assertEquals(new File("/different/../path/path/to/file"), subDir)

    subDir = IOUtils.subDir(new File("/different/./path"), new File("/path/../to/file"))
    Assert.assertEquals(new File("/path/../to/file"), subDir)
  }

  @Test
  def testResetParent = {
    val newFile = IOUtils.resetParent(new File("/new/parent/dir"), new File("/old/parent_dir/file.name"))
    Assert.assertEquals(new File("/new/parent/dir/file.name"), newFile)
  }

  @Test
  def testTempDir = {
    val tempDir = IOUtils.tempDir("Q-Unit-Test")
    Assert.assertTrue(tempDir.exists)
    Assert.assertFalse(tempDir.isFile)
    Assert.assertTrue(tempDir.isDirectory)
    val deleted = tempDir.delete
    Assert.assertTrue(deleted)
    Assert.assertFalse(tempDir.exists)
  }

  @Test
  def testDirLevel = {
    var dir = IOUtils.dirLevel(new File("/path/to/directory"), 1)
    Assert.assertEquals(new File("/path"), dir)

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 2)
    Assert.assertEquals(new File("/path/to"), dir)

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 3)
    Assert.assertEquals(new File("/path/to/directory"), dir)

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 4)
    Assert.assertEquals(new File("/path/to/directory"), dir)
  }

  @Test
  def testAbsolute = {
    var dir = IOUtils.absolute(new File("/path/./to/./directory/."))
    Assert.assertEquals(new File("/path/to/directory"), dir)

    dir = IOUtils.absolute(new File("/"))
    Assert.assertEquals(new File("/"), dir)

    dir = IOUtils.absolute(new File("/."))
    Assert.assertEquals(new File("/"), dir)

    dir = IOUtils.absolute(new File("/././."))
    Assert.assertEquals(new File("/"), dir)

    dir = IOUtils.absolute(new File("/./directory/."))
    Assert.assertEquals(new File("/directory"), dir)

    dir = IOUtils.absolute(new File("/./directory/./"))
    Assert.assertEquals(new File("/directory"), dir)

    dir = IOUtils.absolute(new File("/./directory./"))
    Assert.assertEquals(new File("/directory."), dir)

    dir = IOUtils.absolute(new File("/./.directory/"))
    Assert.assertEquals(new File("/.directory"), dir)
  }
}
