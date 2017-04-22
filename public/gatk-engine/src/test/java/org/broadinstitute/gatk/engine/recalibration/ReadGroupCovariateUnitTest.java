/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.broadinstitute.gatk.engine.recalibration.covariates.ReadGroupCovariate;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
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

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }

    @Test(enabled = true)
    public void testSingleRecord() {
        final String expected = "SAMPLE.1";
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("MY.ID");
        rg.setPlatformUnit(expected);
        runTest(rg, expected, covariate);
    }

    @Test(enabled = true)
    public void testMissingPlatformUnit() {
        final String expected = "MY.7";
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(expected);
        runTest(rg, expected, covariate);
    }

    @Test(enabled = true)
    public void testForceReadgroup() {
        final RecalibrationArgumentCollection forcedRAC = new RecalibrationArgumentCollection();
        forcedRAC.FORCE_READGROUP = "FOO";
        final ReadGroupCovariate forcedCovariate = new ReadGroupCovariate();
        forcedCovariate.initialize(forcedRAC);

        final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("NOT_FOO");
        runTest(rg, "FOO", forcedCovariate);
    }

    private static void runTest(final GATKSAMReadGroupRecord rg, final String expected, final ReadGroupCovariate covariate) {
        GATKSAMRecord read = ReadUtils.createRandomRead(10);
        read.setReadGroup(rg);
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), expected, covariate);

    }

    private static void verifyCovariateArray(final int[][] values, final String expected, final ReadGroupCovariate covariate) {
        for (int[] value : values) {
            String actual = covariate.formatKey(value[0]);
            Assert.assertEquals(actual, expected);
        }
    }

}
