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

package org.broadinstitute.gatk.utils.refdata.tracks;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

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
        File referenceFile = new File(b37KGReference);
        try {
            seq = new CachingIndexedFastaSequenceFile(referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
        genomeLocParser = new GenomeLocParser(seq);

        // We have to disable auto-index creation/locking in the RMDTrackBuilder for tests,
        // as the lock acquisition calls were intermittently hanging on our farm. This unfortunately
        // means that we can't include tests for the auto-index creation feature.
        builder = new RMDTrackBuilder(seq.getSequenceDictionary(),genomeLocParser,null,true,null);
    }

    @Test
    public void testBuilder() {
        Assert.assertTrue(builder.getFeatureManager().getFeatureDescriptors().size() > 0);
    }

    @Test
    public void testDisableAutoIndexGeneration() throws IOException {
        final File unindexedVCF = new File(privateTestDir + "unindexed.vcf");
        final File unindexedVCFIndex = Tribble.indexFile(unindexedVCF);

        Index index = builder.loadIndex(unindexedVCF, new VCFCodec());

        Assert.assertFalse(unindexedVCFIndex.exists());
        Assert.assertNotNull(index);
    }

    @Test
    public void testLoadOnDiskIndex() {
        final File originalVCF = new File(privateTestDir + "vcf4.1.example.vcf");
        final File tempVCFWithCorrectIndex = createTempVCFFileAndIndex(originalVCF, false);
        final File tempVCFIndexFile = Tribble.indexFile(tempVCFWithCorrectIndex);

        final Index index = builder.loadFromDisk(tempVCFWithCorrectIndex, tempVCFIndexFile);

        Assert.assertNotNull(index);
        Assert.assertTrue(tempVCFIndexFile.exists());

        final Index inMemoryIndex = builder.createIndexInMemory(tempVCFWithCorrectIndex, new VCFCodec());
        Assert.assertTrue(index.equalsIgnoreProperties(inMemoryIndex));
    }

    @Test
    public void testLoadOnDiskOutdatedIndex() {
        final File originalVCF = new File(privateTestDir + "vcf4.1.example.vcf");
        final File tempVCFWithOutdatedIndex = createTempVCFFileAndIndex(originalVCF, true);
        final File tempVCFIndexFile = Tribble.indexFile(tempVCFWithOutdatedIndex);

        final Index index = builder.loadFromDisk(tempVCFWithOutdatedIndex, tempVCFIndexFile);

        // loadFromDisk() should return null to indicate that the index is outdated and should not be used,
        // but should not delete the index since our builder has disableAutoIndexCreation set to true
        Assert.assertNull(index);
        Assert.assertTrue(tempVCFIndexFile.exists());
    }

    /**
     * Create a temporary vcf file and an associated index file, which may be set to be out-of-date
     * relative to the vcf
     *
     * @param vcfFile the vcf file
     * @param createOutOfDateIndex if true, ensure that the temporary vcf file is modified after the index
     * @return a file pointing to the new tmp location, with accompanying index
     */
    private File createTempVCFFileAndIndex( final File vcfFile, final boolean createOutOfDateIndex ) {
        try {
            final File tmpFile = createTempFile("RMDTrackBuilderUnitTest", "");
            final File tmpIndex = Tribble.indexFile(tmpFile);
            tmpIndex.deleteOnExit();

            copyFile(vcfFile, tmpFile);
            final Index inMemoryIndex = builder.createIndexInMemory(tmpFile, new VCFCodec());
            final LittleEndianOutputStream indexOutputStream = new LittleEndianOutputStream(new FileOutputStream(tmpIndex));

            // If requested, modify the tribble file after the index. Otherwise, modify the index last.
            if ( createOutOfDateIndex ) {
                inMemoryIndex.write(indexOutputStream);
                indexOutputStream.close();
                Thread.sleep(2000);
                copyFile(vcfFile, tmpFile);
            }
            else {
                copyFile(vcfFile, tmpFile);
                Thread.sleep(2000);
                inMemoryIndex.write(indexOutputStream);
                indexOutputStream.close();
            }

            return tmpFile;
        } catch (IOException e) {
            Assert.fail("Unable to create temperary file");
        } catch (InterruptedException e) {
            Assert.fail("Somehow our thread got interrupted");
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

