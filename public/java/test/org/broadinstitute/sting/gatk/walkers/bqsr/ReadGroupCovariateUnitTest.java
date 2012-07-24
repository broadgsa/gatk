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
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), expected);

    }

    private void verifyCovariateArray(int[][] values, String expected) {
        for (int[] value : values) {
            String actual = covariate.formatKey(value[0]);
            Assert.assertEquals(actual, expected);
        }
    }

}
