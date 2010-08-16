/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks.builders;


import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broad.tribble.index.Index;
import org.broad.tribble.vcf.VCFCodec;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.nio.channels.FileChannel;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class TribbleRMDTrackBuilderUnitTest
 *         <p/>
 *         Testing out the builder for tribble Tracks
 */
public class TribbleRMDTrackBuilderUnitTest extends BaseTest {
    private TribbleRMDTrackBuilder builder;

    @Before
    public void setup() {
        builder = new TribbleRMDTrackBuilder();
        IndexedFastaSequenceFile seq = new IndexedFastaSequenceFile(new File(b36KGReference));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Test
    public void testBuilder() {
        Map<String, Class> classes = builder.getAvailableTrackNamesAndTypes();
        Assert.assertTrue(classes.size() > 0);
    }

    @Test
    // in this test, the index exists, but is out of date.
    public void testBuilderIndexUnwriteable() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/read_only/relic.vcf");
        try {
            builder.loadIndex(vcfFile, new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we didn't write the file (check that it's timestamp is within bounds)
        //System.err.println(new File(vcfFile + TribbleRMDTrackBuilder.indexExtension).lastModified());
        Assert.assertTrue(Math.abs(1279591752000l - new File(vcfFile + TribbleRMDTrackBuilder.indexExtension).lastModified()) < 100);

    }

    // we have a good index file, in a read-only dir. This would cause the previous version to remake the index; make
    // sure we don't do this
    @Test
    public void testDirIsLockedIndexFromDisk() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/read_only/good_index.vcf");
        File vcfFileIndex = new File(validationDataLocation + "/ROD_validation/read_only/good_index.vcf.idx");
        Index ind = null;
        try {
            ind = builder.attemptIndexFromDisk(vcfFile,new VCFCodec(),vcfFileIndex,new FSLockWithShared(vcfFile));
        } catch (IOException e) {
            Assert.fail("We weren't expecting an exception -> " + e.getMessage());
        }
        // make sure we get back a null index; i.e. we can't load the index from disk
        Assert.assertTrue(ind == null);
    }



    @Test
    public void testBuilderIndexDirectoryUnwritable() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/read_only/no_index.vcf");
        File vcfFileIndex = new File(validationDataLocation + "/ROD_validation/read_only/no_index.vcf.idx");

        Index ind = null;
        try {
            ind = builder.loadIndex(vcfFile, new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we didn't write the file (check that it's timestamp is within bounds)
        Assert.assertTrue(!vcfFileIndex.exists());
        Assert.assertTrue(ind != null);

    }


    @Test
    public void testGenerateIndexForUnindexedFile() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/always_reindex.vcf");
        File vcfFileIndex = new File(validationDataLocation + "/ROD_validation/always_reindex.vcf.idx");

        // if we can't write to the directory, don't fault the tester, just pass
        if (!vcfFileIndex.getParentFile().canWrite()) {
            logger.warn("Unable to run test testGenerateIndexForUnindexedFile: unable to write to dir " + vcfFileIndex.getParentFile());
            return;
        }
        // clean-up our test, and previous tests that may have written the file
        vcfFileIndex.deleteOnExit();
        if (vcfFileIndex.exists()) vcfFileIndex.delete();

        try {
            builder.loadIndex(vcfFile, new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we wrote the file
        Assert.assertTrue(vcfFileIndex.exists());
    }


    // test to make sure we delete the index and regenerate if it's out of date
    //@Test
    public void testBuilderIndexOutOfDate() {
        File vcfFile = createOutofDateIndexFile(new File(validationDataLocation + "/ROD_validation/newerTribbleTrack.vcf"));
        try {
            builder.loadIndex(vcfFile, new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        //System.err.println("index : " + new File(vcfFile + ".idx").lastModified());
        // System.err.println("vcf : " + vcfFile.lastModified());

        // make sure that we removed and updated the index
        Assert.assertTrue("VCF file index wasn't updated", new File(vcfFile + ".idx").lastModified() > vcfFile.lastModified());
    }

    // test to make sure we delete the index and regenerate if it's out of date
    // @Test
    public void testBuilderIndexGoodDate() {
        File vcfFile = createCorrectDateIndexFile(new File(validationDataLocation + "/ROD_validation/newerTribbleTrack.vcf"));
        Long indexTimeStamp = new File(vcfFile.getAbsolutePath() + ".idx").lastModified();
        try {
            builder.loadIndex(vcfFile, new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        //System.err.println("index : " + new File(vcfFile + ".idx").lastModified());
        //System.err.println("old : " + indexTimeStamp);

        // make sure that we removed and updated the index
        Assert.assertTrue("Fail: index file was modified", new File(vcfFile + ".idx").lastModified() == indexTimeStamp);
    }

    /**
     * create a temporary file and an associated out of date index file
     *
     * @param tribbleFile the tribble file
     * @return a file pointing to the new tmp location, with out of date index
     */
    private File createOutofDateIndexFile(File tribbleFile) {
        try {
            // first copy the tribble file to a temperary file
            File tmpFile = File.createTempFile("TribbleUnitTestFile", "");
            tmpFile.deleteOnExit();

            logger.info("creating temp file " + tmpFile);
            // create a fake index, before we copy so it's out of date
            File tmpIndex = new File(tmpFile.getAbsolutePath() + ".idx");
            tmpIndex.deleteOnExit();

            // sleep, to make sure the timestamps are different
            Thread.sleep(2000);

            // copy the vcf (tribble) file to the tmp file location
            copyFile(tribbleFile, tmpFile);

            // sleep again, to make sure the timestamps are different (vcf vrs updated index file)
            Thread.sleep(2000);

            return tmpFile;

        } catch (IOException e) {
            Assert.fail("Fail: Unable to create temperary file");
        } catch (InterruptedException e) {
            Assert.fail("Fail: Somehow our thread got interupted");
        }
        return null;
    }

    /**
     * create a temporary file and an associated out of date index file
     *
     * @param tribbleFile the tribble file
     * @return a file pointing to the new tmp location, with out of date index
     */
    private File createCorrectDateIndexFile(File tribbleFile) {
        try {
            // first copy the tribble file to a temperary file
            File tmpFile = File.createTempFile("TribbleUnitTestFile", "");
            tmpFile.deleteOnExit();
            logger.info("creating temp file " + tmpFile);

            // copy the vcf (tribble) file to the tmp file location
            copyFile(tribbleFile, tmpFile);

            // sleep again, to make sure the timestamps are different (vcf vrs updated index file)
            Thread.sleep(2000);

            // create a fake index, before we copy so it's out of date
            File tmpIndex = new File(tmpFile.getAbsolutePath() + ".idx");
            tmpIndex.deleteOnExit();

            // copy the vcf (tribble) file to the tmp file location
            copyFile(new File(tribbleFile + ".idx"), tmpIndex);

            return tmpFile;

        } catch (IOException e) {
            Assert.fail("Unable to create temperary file");
        } catch (InterruptedException e) {
            Assert.fail("Somehow our thread got interupted");
        }
        return null;
    }

    /**
     * copy a file, from http://www.exampledepot.com/egs/java.nio/File2File.html
     *
     * @param srFile the source file
     * @param dtFile the destination file
     */
    private static void copyFile(File srFile, File dtFile) {
        try {
            // Create channel on the source
            FileChannel srcChannel = new FileInputStream(srFile).getChannel();

            // Create channel on the destination
            FileChannel dstChannel = new FileOutputStream(dtFile).getChannel();

            // Copy file contents from source to destination
            dstChannel.transferFrom(srcChannel, 0, srcChannel.size());

            // Close the channels
            srcChannel.close();
            dstChannel.close();
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("Unable to process copy " + e.getMessage());
        }
    }

}

