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

package org.broadinstitute.gatk.engine;

import htsjdk.samtools.SAMFileReader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Eric Banks
 * @since 7/18/12
 */
public class CommandLineGATKUnitTest extends BaseTest {

    @Test(enabled = true)
    public void testSamTextFileError1() {
        final File samFile = new File(publicTestDir + "testfile.sam");
        final File indexFile = new File(publicTestDir + "HiSeq.1mb.1RG.bai");
        try {
            final SAMFileReader reader = new SAMFileReader(samFile, indexFile, false);

            // we shouldn't get here
            Assert.fail("We should have exceptioned out when trying to create a reader with an index for a textual SAM file");
        } catch (RuntimeException e) {
            Assert.assertTrue(e.getMessage().indexOf(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_1) != -1);
        }
    }

    @Test(enabled = true)
    public void testSamTextFileError2() {
        File samFile = new File(publicTestDir + "testfile.sam");
        try {
            final SAMFileReader reader = new SAMFileReader(samFile);
            reader.getFilePointerSpanningReads();

            // we shouldn't get here
            Assert.fail("We should have exceptioned out when trying to call getFilePointerSpanningReads() for a textual SAM file");
        } catch (RuntimeException e) {
            Assert.assertTrue(e.getMessage().indexOf(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_2) != -1);
        }
    }
}
