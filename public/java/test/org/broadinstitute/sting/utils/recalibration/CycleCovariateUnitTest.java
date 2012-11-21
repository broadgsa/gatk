package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.covariates.CycleCovariate;
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

    @Test(enabled = true)
    public void testSimpleCycles() {
        short readLength = 10;
        GATKSAMRecord read = ReadUtils.createRandomRead(readLength);
        read.setReadPairedFlag(true);
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), 1, (short) 1);

        read.setReadNegativeStrandFlag(true);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), readLength, -1);

        read.setSecondOfPairFlag(true);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), -readLength, 1);

        read.setReadNegativeStrandFlag(false);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), -1, -1);
    }

    private void verifyCovariateArray(int[][] values, int init, int increment) {
        for (short i = 0; i < values.length; i++) {
            short actual = Short.decode(covariate.formatKey(values[i][0]));
            int expected = init + (increment * i);
            Assert.assertEquals(actual, expected);
        }
    }

    @Test(enabled = true, expectedExceptions={UserException.class})
    public void testMoreThanMaxCycleFails() {
        int readLength = RAC.MAXIMUM_CYCLE_VALUE + 1;
        GATKSAMRecord read = ReadUtils.createRandomRead(readLength);
        read.setReadPairedFlag(true);
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
    }

    @Test(enabled = true)
    public void testMaxCyclePasses() {
        int readLength = RAC.MAXIMUM_CYCLE_VALUE;
        GATKSAMRecord read = ReadUtils.createRandomRead(readLength);
        read.setReadPairedFlag(true);
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
    }
}
