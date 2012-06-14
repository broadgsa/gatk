package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * @author Mauricio Carneiro
 * @since 3/1/12
 */
public class ReadGroupCovariateUnitTest {
    ReadGroupCovariate covariate;
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ReadGroupCovariate();
        covariate.initialize(RAC);
    }

    @Test(enabled = true)
    public void testSingleRecord() {
        final String expected = "SAMPLE.1";
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("MY.ID");
        rg.setPlatformUnit(expected);
        runTest(rg, expected);
    }

    @Test(enabled = true)
    public void testMissingPlatformUnit() {
        final String expected = "MY.7";
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(expected);
        runTest(rg, expected);
    }

    private void runTest(GATKSAMReadGroupRecord rg, String expected) {
        GATKSAMRecord read = ReadUtils.createRandomRead(10);
        read.setReadGroup(rg);
        CovariateValues values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), expected);

    }

    private void verifyCovariateArray(long[] values, String expected) {
        for (Long value : values) {
            String actual = covariate.formatKey(value);
            Assert.assertEquals(actual, expected);
        }
    }

}
