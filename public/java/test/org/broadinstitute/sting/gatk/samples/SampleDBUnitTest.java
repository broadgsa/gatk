package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;

import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleDBUnitTest extends BaseTest {
    // this empty header used to instantiate sampledatasource objects
    private static SAMFileHeader header = new SAMFileHeader();

    // all the test sample files are located here
    private String sampleFilesDir = validationDataLocation +  "samples/";

    // make sure samples are created from the SAM file correctly
    @Test()
    public void loadSAMSamplesTest() {
        //SampleDB s = new SampleDB(header);
    }
}
