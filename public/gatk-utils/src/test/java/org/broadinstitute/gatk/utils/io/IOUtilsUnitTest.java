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

package org.broadinstitute.gatk.utils.io;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
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
            IOUtils.writeResource(new Resource("testProperties.properties", null), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testWriteSystemTempFile() throws IOException {
        File temp = IOUtils.writeTempResource(new Resource("testProperties.properties", null));
        try {
            Assert.assertTrue(temp.getName().startsWith("testProperties"), "File does not start with 'testProperties.': " + temp);
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
            IOUtils.writeResource(new Resource("/testProperties.properties", IOUtils.class), temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test
    public void testWriteRelativeTempFile() throws IOException {
        File temp = IOUtils.writeTempResource(new Resource("/testProperties.properties", IOUtils.class));
        try {
            Assert.assertTrue(temp.getName().startsWith("testProperties"), "File does not start with 'testProperties.': " + temp);
            Assert.assertTrue(temp.getName().endsWith(".properties"), "File does not end with '.properties': " + temp);
        } finally {
            FileUtils.deleteQuietly(temp);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingRelativeFile() throws IOException {
        File temp = createTempFile("temp.", ".properties");
        try {
            // Looking for /org/broadinstitute/gatk/utils/file/GATKText.properties
            IOUtils.writeResource(new Resource("GATKText.properties", IOUtils.class), temp);
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

    @Test
    public void testIsSpecialFile() {
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/null")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/full")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stdout")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stderr")));
        Assert.assertFalse(IOUtils.isSpecialFile(null));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/home/user/my.file")));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/devfake/null")));
    }

    @DataProvider( name = "ByteArrayIOTestData")
    public Object[][] byteArrayIOTestDataProvider() {
        return new Object[][] {
            // file size, read buffer size
            { 0,     4096 },
            { 1,     4096 },
            { 2000,  4096 },
            { 4095,  4096 },
            { 4096,  4096 },
            { 4097,  4096 },
            { 6000,  4096 },
            { 8191,  4096 },
            { 8192,  4096 },
            { 8193,  4096 },
            { 10000, 4096 }
        };
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadFileIntoByteArray ( int fileSize, int readBufferSize ) throws Exception {
        File tempFile = createTempFile(String.format("testWriteThenReadFileIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToFile(dataWritten, tempFile);
        byte[] dataRead = IOUtils.readFileIntoByteArray(tempFile, readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertTrue(Arrays.equals(dataRead, dataWritten));
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadStreamIntoByteArray ( int fileSize, int readBufferSize ) throws Exception {
        File tempFile = createTempFile(String.format("testWriteThenReadStreamIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToStream(dataWritten, new FileOutputStream(tempFile));
        byte[] dataRead = IOUtils.readStreamIntoByteArray(new FileInputStream(tempFile), readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertTrue(Arrays.equals(dataRead, dataWritten));
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testReadNonExistentFileIntoByteArray() {
        File nonExistentFile = new File("djfhsdkjghdfk");
        Assert.assertFalse(nonExistentFile.exists());

        IOUtils.readFileIntoByteArray(nonExistentFile);
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testReadNullStreamIntoByteArray() {
        IOUtils.readStreamIntoByteArray(null);
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testReadStreamIntoByteArrayInvalidBufferSize() throws Exception {
        IOUtils.readStreamIntoByteArray(new FileInputStream(createTempFile("testReadStreamIntoByteArrayInvalidBufferSize", "tmp")),
                                        -1);
    }

    @Test( expectedExceptions = UserException.CouldNotCreateOutputFile.class )
    public void testWriteByteArrayToUncreatableFile() {
        IOUtils.writeByteArrayToFile(new byte[]{0}, new File("/dev/foo/bar"));
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testWriteNullByteArrayToFile() {
        IOUtils.writeByteArrayToFile(null, createTempFile("testWriteNullByteArrayToFile", "tmp"));
    }

    @Test( expectedExceptions = ReviewedGATKException.class )
    public void testWriteByteArrayToNullStream() {
        IOUtils.writeByteArrayToStream(new byte[]{0}, null);
    }

    private byte[] getDeterministicRandomData ( int size ) {
        Utils.resetRandomGenerator();
        Random rand = Utils.getRandomGenerator();

        byte[] randomData = new byte[size];
        rand.nextBytes(randomData);

        return randomData;
    }
}
