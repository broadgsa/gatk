package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.StingException;

import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleDataSourceUnitTest extends BaseTest {
    // this empty header used to instantiate sampledatasource objects
    private static SAMFileHeader header = new SAMFileHeader();

    // all the test sample files are located here
    private String sampleFilesDir = validationDataLocation +  "samples/";

    // make sure samples are created from the SAM file correctly
    @Test()
    public void loadSAMSamplesTest() {
        SampleDataSource s = new SampleDataSource(header, null);
    }
}
