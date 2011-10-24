package org.broadinstitute.sting.utils.io;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.sting.BaseTest;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

public class IOUtilsUnitTest extends BaseTest {
  @Test
  public void testGoodTempDir() {
    IOUtils.checkTempDir(new File("/tmp/queue"));
  }

  @Test(expectedExceptions=UserException.BadTmpDir.class)
  public void testBadTempDir() {
    IOUtils.checkTempDir(new File("/tmp"));
  }

  @Test
  public void testAbsoluteSubDir() {
    File subDir = IOUtils.absolute(new File("."), new File("/path/to/file"));
    Assert.assertEquals(subDir, new File("/path/to/file"));

    subDir = IOUtils.absolute(new File("/different/path"), new File("/path/to/file"));
    Assert.assertEquals(subDir, new File("/path/to/file"));

    subDir = IOUtils.absolute(new File("/different/path"), new File("."));
    Assert.assertEquals(subDir, new File("/different/path"));
  }

  @Test
  public void testRelativeSubDir() throws IOException {
    File subDir = IOUtils.absolute(new File("."), new File("path/to/file"));
    Assert.assertEquals(subDir.getCanonicalFile(), new File("path/to/file").getCanonicalFile());

    subDir = IOUtils.absolute(new File("/different/path"), new File("path/to/file"));
    Assert.assertEquals(subDir, new File("/different/path/path/to/file"));
  }

  @Test
  public void testDottedSubDir() throws IOException {
    File subDir = IOUtils.absolute(new File("."), new File("path/../to/file"));
    Assert.assertEquals(subDir.getCanonicalFile(), new File("path/../to/./file").getCanonicalFile());

    subDir = IOUtils.absolute(new File("."), new File("/path/../to/file"));
    Assert.assertEquals(subDir, new File("/path/../to/file"));

    subDir = IOUtils.absolute(new File("/different/../path"), new File("path/to/file"));
    Assert.assertEquals(subDir, new File("/different/../path/path/to/file"));

    subDir = IOUtils.absolute(new File("/different/./path"), new File("/path/../to/file"));
    Assert.assertEquals(subDir, new File("/path/../to/file"));
  }

  @Test
  public void testTempDir() {
    File tempDir = IOUtils.tempDir("Q-Unit-Test", "", new File("queueTempDirToDelete"));
    Assert.assertTrue(tempDir.exists());
    Assert.assertFalse(tempDir.isFile());
    Assert.assertTrue(tempDir.isDirectory());
    boolean deleted = IOUtils.tryDelete(tempDir);
    Assert.assertTrue(deleted);
    Assert.assertFalse(tempDir.exists());
  }

  @Test
  public void testDirLevel() {
    File dir = IOUtils.dirLevel(new File("/path/to/directory"), 1);
    Assert.assertEquals(dir, new File("/path"));

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 2);
    Assert.assertEquals(dir, new File("/path/to"));

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 3);
    Assert.assertEquals(dir, new File("/path/to/directory"));

    dir = IOUtils.dirLevel(new File("/path/to/directory"), 4);
    Assert.assertEquals(dir, new File("/path/to/directory"));
  }

  @Test
  public void testAbsolute() {
    File dir = IOUtils.absolute(new File("/path/./to/./directory/."));
    Assert.assertEquals(dir, new File("/path/to/directory"));

    dir = IOUtils.absolute(new File("/"));
    Assert.assertEquals(dir, new File("/"));

    dir = IOUtils.absolute(new File("/."));
    Assert.assertEquals(dir, new File("/"));

    dir = IOUtils.absolute(new File("/././."));
    Assert.assertEquals(dir, new File("/"));

    dir = IOUtils.absolute(new File("/./directory/."));
    Assert.assertEquals(dir, new File("/directory"));

    dir = IOUtils.absolute(new File("/./directory/./"));
    Assert.assertEquals(dir, new File("/directory"));

    dir = IOUtils.absolute(new File("/./directory./"));
    Assert.assertEquals(dir, new File("/directory."));

    dir = IOUtils.absolute(new File("/./.directory/"));
    Assert.assertEquals(dir, new File("/.directory"));
  }

  @Test
  public void testTail() throws IOException {
    List<String> lines = Arrays.asList(
            "chr18_random	4262	3154410390	50	51",
            "chr19_random	301858	3154414752	50	51",
            "chr21_random	1679693	3154722662	50	51",
            "chr22_random	257318	3156435963	50	51",
            "chrX_random	1719168	3156698441	50	51");
    List<String> tail = IOUtils.tail(new File(BaseTest.hg18Reference + ".fai"), 5);
    Assert.assertEquals(tail.size(), 5);
    for (int i = 0; i < 5; i++)
      Assert.assertEquals(tail.get(i), lines.get(i));
  }

    @Test
    public void testWriteSystemFile() throws IOException {
        File temp = createTempFile("temp.", ".properties");
        try {
            IOUtils.writeResource(new Resource("StingText.properties", null), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testWriteSystemTempFile() throws IOException {
        File temp = IOUtils.writeTempResource(new Resource("StingText.properties", null));
        try {
            Assert.assertTrue(temp.getName().startsWith("StingText"), "File does not start with 'StingText.': " + temp);
            Assert.assertTrue(temp.getName().endsWith(".properties"), "File does not end with '.properties': " + temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingSystemFile() throws IOException {
        File temp = createTempFile("temp.", ".properties");
        try {
            IOUtils.writeResource(new Resource("MissingStingText.properties", null), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testWriteRelativeFile() throws IOException {
        File temp = createTempFile("temp.", ".properties");
        try {
            IOUtils.writeResource(new Resource("/StingText.properties", IOUtils.class), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testWriteRelativeTempFile() throws IOException {
        File temp = IOUtils.writeTempResource(new Resource("/StingText.properties", IOUtils.class));
        try {
            Assert.assertTrue(temp.getName().startsWith("StingText"), "File does not start with 'StingText.': " + temp);
            Assert.assertTrue(temp.getName().endsWith(".properties"), "File does not end with '.properties': " + temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingRelativeFile() throws IOException {
        File temp = createTempFile("temp.", ".properties");
        try {
            // Looking for /org/broadinstitute/sting/utils/file/StingText.properties
            IOUtils.writeResource(new Resource("StingText.properties", IOUtils.class), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testResourceProperties() {
        Resource resource = new Resource("foo", Resource.class);
        Assert.assertEquals(resource.getPath(), "foo");
        Assert.assertEquals(resource.getRelativeClass(), Resource.class);
    }
}
