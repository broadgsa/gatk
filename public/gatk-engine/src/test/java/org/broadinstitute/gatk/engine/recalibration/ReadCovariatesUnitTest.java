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

import org.broadinstitute.gatk.engine.recalibration.covariates.*;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.recalibration.EventType;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Random;

/**
 * @author carneiro
 * @since 4/21/12
 */
public class ReadCovariatesUnitTest {

    @BeforeMethod
    public void init() {
        ReadCovariates.clearKeysCache();
    }

    @Test(enabled = false)
    public void testCovariateGeneration() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        final String RGID = "id";

        ReadGroupCovariate rgCov = new ReadGroupCovariate();
        QualityScoreCovariate qsCov = new QualityScoreCovariate();
        ContextCovariate coCov = new ContextCovariate();
        CycleCovariate cyCov = new CycleCovariate();

        rgCov.initialize(RAC);
        qsCov.initialize(RAC);
        coCov.initialize(RAC);
        cyCov.initialize(RAC);

        Covariate[] requestedCovariates = new Covariate[4];
        requestedCovariates[0] = rgCov;
        requestedCovariates[1] = qsCov;
        requestedCovariates[2] = coCov;
        requestedCovariates[3] = cyCov;

        final int NUM_READS = 100;
        final Random rnd = Utils.getRandomGenerator();

        final String[] readGroups = {"RG1", "RG2", "RGbla"};
        for (int idx = 0; idx < NUM_READS; idx++) {
            for (final String rgs : readGroups) {
                final int length = 10 + rnd.nextInt(100); // random read length, at least 10 bp long
                final GATKSAMRecord read = ReadUtils.createRandomRead(length, false);
                final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(rgs);
                rg.setPlatform("illumina");
                read.setReadGroup(rg);
                read.setReadNegativeStrandFlag(rnd.nextBoolean());
                final byte[] mQuals = read.getBaseQualities(EventType.BASE_SUBSTITUTION);
                final byte[] iQuals = read.getBaseQualities(EventType.BASE_INSERTION);
                final byte[] dQuals = read.getBaseQualities(EventType.BASE_DELETION);
                ReadCovariates rc = RecalUtils.computeCovariates(read, requestedCovariates);

                // check that the length is correct
                Assert.assertEquals(rc.getMismatchesKeySet().length, length);
                Assert.assertEquals(rc.getInsertionsKeySet().length, length);
                Assert.assertEquals(rc.getDeletionsKeySet().length, length);

                for (int i = 0; i < length; i++) {
                    // check that read group is always the same
                    Assert.assertEquals(rgCov.formatKey(rc.getMismatchesKeySet(i)[0]), rgs);
                    Assert.assertEquals(rgCov.formatKey(rc.getInsertionsKeySet(i)[0]), rgs);
                    Assert.assertEquals(rgCov.formatKey(rc.getDeletionsKeySet(i)[0]),  rgs);

                    // check quality score
                    Assert.assertEquals(qsCov.formatKey(rc.getMismatchesKeySet(i)[1]), "" + mQuals[i]);
                    Assert.assertEquals(qsCov.formatKey(rc.getInsertionsKeySet(i)[1]), "" + iQuals[i]);
                    Assert.assertEquals(qsCov.formatKey(rc.getDeletionsKeySet(i)[1]),  "" + dQuals[i]);

                    // check context
                    Assert.assertEquals(coCov.formatKey(rc.getMismatchesKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.MISMATCHES_CONTEXT_SIZE));
                    Assert.assertEquals(coCov.formatKey(rc.getInsertionsKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE));
                    Assert.assertEquals(coCov.formatKey(rc.getDeletionsKeySet(i)[2]),  ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE));

                    // check cycle
                    Assert.assertEquals(cyCov.formatKey(rc.getMismatchesKeySet(i)[3]), "" + (i+1));
                    Assert.assertEquals(cyCov.formatKey(rc.getInsertionsKeySet(i)[3]), "" + (i+1));
                    Assert.assertEquals(cyCov.formatKey(rc.getDeletionsKeySet(i)[3]),  "" + (i+1));
                }

            }

        }

    }

}
