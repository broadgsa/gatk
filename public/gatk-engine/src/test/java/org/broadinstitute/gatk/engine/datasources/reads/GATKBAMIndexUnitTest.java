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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Test basic functionality in the GATK's implementation of the BAM index classes.
 */
public class GATKBAMIndexUnitTest extends BaseTest {
    private static File bamFile = new File(validationDataLocation+"MV1994.selected.bam");

    /**
     * Index file forming the source of all unit tests.
     */
    private static File bamIndexFile = new File(validationDataLocation+"MV1994.selected.bam.bai");

    /**
     * Storage for the index itself.
     */
    private GATKBAMIndex bamIndex;

    /**
     * Sequences.
     */
    private SAMSequenceDictionary sequenceDictionary;


    @BeforeClass
    public void init() throws FileNotFoundException {
        SAMFileReader reader = new SAMFileReader(bamFile);
        this.sequenceDictionary = reader.getFileHeader().getSequenceDictionary();
        reader.close();
        
        bamIndex = new GATKBAMIndex(bamIndexFile, sequenceDictionary);
    }

    @Test
    public void testNumberAndSizeOfIndexLevels() {
        // The correct values for this test are pulled directly from the
        // SAM Format Specification v1.3-r882, Section 4.1.1, last paragraph.
        Assert.assertEquals(GATKBAMIndex.getNumIndexLevels(),6,"Incorrect number of levels in BAM index");

        // Level 0
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(0),0);
        Assert.assertEquals(bamIndex.getLevelSize(0),1);

        // Level 1
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(1),1);
        Assert.assertEquals(bamIndex.getLevelSize(1),8-1+1);

        // Level 2
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(2),9);
        Assert.assertEquals(bamIndex.getLevelSize(2),72-9+1);

        // Level 3
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(3),73);
        Assert.assertEquals(bamIndex.getLevelSize(3),584-73+1);

        // Level 4
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(4),585);
        Assert.assertEquals(bamIndex.getLevelSize(4),4680-585+1);

        // Level 5                                
        Assert.assertEquals(GATKBAMIndex.getFirstBinInLevel(5),4681);
        Assert.assertEquals(bamIndex.getLevelSize(5),37448-4681+1);
    }

    @Test( expectedExceptions = UserException.MalformedFile.class )
    public void testDetectTruncatedBamIndexWordBoundary() {
        GATKBAMIndex index = new GATKBAMIndex(new File(privateTestDir + "truncated_at_word_boundary.bai"), sequenceDictionary);
        index.readReferenceSequence(0);
    }

    @Test( expectedExceptions = UserException.MalformedFile.class )
    public void testDetectTruncatedBamIndexNonWordBoundary() {
        GATKBAMIndex index = new GATKBAMIndex(new File(privateTestDir + "truncated_at_non_word_boundary.bai"), sequenceDictionary);
        index.readReferenceSequence(0);
    }

}
