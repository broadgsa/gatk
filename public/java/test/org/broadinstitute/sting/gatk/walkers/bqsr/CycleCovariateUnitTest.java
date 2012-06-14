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
public class CycleCovariateUnitTest {
    CycleCovariate covariate;
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new CycleCovariate();
        covariate.initialize(RAC);
    }

    @Test(enabled = false)
    public void testSimpleCycles() {
        short readLength = 10;
        GATKSAMRecord read = ReadUtils.createRandomRead(readLength);
        read.setReadPairedFlag(true);
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        CovariateValues values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), 1, (short) 1);

        read.setReadNegativeStrandFlag(true);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), readLength, -1);

        read.setSecondOfPairFlag(true);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), -readLength, 1);

        read.setReadNegativeStrandFlag(false);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), -1, -1);
    }

    private void verifyCovariateArray(long[] values, int init, int increment) {
        for (short i = 0; i < values.length; i++) {
            short actual = Short.decode(covariate.formatKey(values[i]));
            int expected = init + (increment * i);
            //            System.out.println(String.format("%d: %d, %d", i, actual, expected));
            Assert.assertEquals(actual, expected);
        }
    }

}
