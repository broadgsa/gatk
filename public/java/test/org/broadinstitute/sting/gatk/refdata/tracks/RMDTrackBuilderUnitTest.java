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

package org.broadinstitute.sting.gatk.refdata.tracks;


import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.file.FSLockWithShared;

import org.testng.annotations.BeforeMethod;

import org.testng.annotations.Test;

import java.io.*;
import java.nio.channels.FileChannel;


/**
 * @author aaron
 *         <p/>
 *         Class RMDTrackBuilderUnitTest
 *         <p/>
 *         Testing out the builder for tribble Tracks
 */
public class RMDTrackBuilderUnitTest extends BaseTest {
    private RMDTrackBuilder builder;
    private IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeMethod
    public void setup() {
        File referenceFile = new File(b36KGReference);
        try {
            seq = new CachingIndexedFastaSequenceFile(referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
        genomeLocParser = new GenomeLocParser(seq);
        builder = new RMDTrackBuilder(seq.getSequenceDictionary(),genomeLocParser,null);
    }

    @Test
    public void testBuilder() {
        Assert.assertTrue(builder.getFeatureManager().getFeatureDescriptors().size() > 0);
    }

    @Test
    // in this test, the index exists, but is out of date.
    public void testBuilderIndexUnwriteable() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/read_only/relic.vcf");
        try {
            builder.loadIndex(vcfFile, new VCF3Codec());
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we didn't write the file (check that it's timestamp is within bounds)
        //System.err.println(new File(vcfFile + RMDTrackBuilder.indexExtension).lastModified());
        Assert.assertTrue(Math.abs(1279591752000l - Tribble.indexFile(vcfFile).lastModified()) < 100);

    }

    // we have a good index file, in a read-only dir. This would cause the previous version to remake the index; make
    // sure we don't do this
    @Test
    public void testDirIsLockedIndexFromDisk() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/read_only/good_index.vcf");
        File vcfFileIndex = Tribble.indexFile(vcfFile);
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
        File vcfFileIndex = Tribble.indexFile(vcfFile);

        Index ind = null;
        try {
            ind = builder.loadIndex(vcfFile, new VCF3Codec());
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
        File vcfFileIndex = Tribble.indexFile(vcfFile);

        // if we can't write to the directory, don't fault the tester, just pass
        if (!vcfFileIndex.getParentFile().canWrite()) {
            logger.warn("Unable to run test testGenerateIndexForUnindexedFile: unable to write to dir " + vcfFileIndex.getParentFile());
            return;
        }
        // clean-up our test, and previous tests that may have written the file
        vcfFileIndex.deleteOnExit();
        if (vcfFileIndex.exists()) vcfFileIndex.delete();

        try {
            builder.loadIndex(vcfFile, new VCFCodec());
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we wrote the file
        Assert.assertTrue(vcfFileIndex.exists());
    }


    // test to make sure we get a full sequence dictionary from the VCF (when we set the dictionary in the builder)
    @Test
    public void testBuilderIndexSequenceDictionary() {
        File vcfFile = createCorrectDateIndexFile(new File(validationDataLocation + "/ROD_validation/newerTribbleTrack.vcf"));
        Long indexTimeStamp = Tribble.indexFile(vcfFile).lastModified();
        try {
            Index idx = builder.loadIndex(vcfFile, new VCFCodec());
            // catch any exception; this call should pass correctly
            SAMSequenceDictionary dict =  IndexDictionaryUtils.getSequenceDictionaryFromProperties(idx);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }

        // make sure that we removed and updated the index
        Assert.assertTrue(Tribble.indexFile(vcfFile).lastModified() >= indexTimeStamp,"Fail: index file was modified");
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
            File tmpIndex = Tribble.indexFile(tmpFile);
            tmpIndex.deleteOnExit();

            // copy the vcf (tribble) file to the tmp file location
            copyFile(Tribble.indexFile(tribbleFile), tmpIndex);

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

