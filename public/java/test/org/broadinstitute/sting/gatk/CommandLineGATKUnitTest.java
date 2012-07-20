package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.BaseTest;
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
