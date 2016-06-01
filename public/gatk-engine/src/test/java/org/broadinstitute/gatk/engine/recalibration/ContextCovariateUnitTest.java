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

import org.broadinstitute.gatk.engine.recalibration.covariates.ContextCovariate;
import org.broadinstitute.gatk.engine.recalibration.covariates.Covariate;
import org.broadinstitute.gatk.utils.clipping.ClippingRepresentation;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
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
public class ContextCovariateUnitTest {
    ContextCovariate covariate;
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ContextCovariate();
        covariate.initialize(RAC);
    }

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }

    @Test(enabled = true)
    public void testSimpleContexts() {
        GATKSAMRecord read = ReadUtils.createRandomRead(1000);
        GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);

        verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getDeletionsKeySet(),  RAC.INDELS_CONTEXT_SIZE,  clippedRead, covariate);
    }

    public static void verifyCovariateArray(int[][] values, int contextSize, GATKSAMRecord read, Covariate contextCovariate) {
        for (int i = 0; i < values.length; i++)
            Assert.assertEquals(contextCovariate.formatKey(values[i][0]), expectedContext(read, i, contextSize));

    }

    public static String expectedContext (GATKSAMRecord read, int offset, int contextSize) {
        final String bases = stringFrom(read.getReadBases());
        String expectedContext = null;
        if (offset - contextSize + 1 >= 0) {
            String context = bases.substring(offset - contextSize + 1, offset + 1);
            if (!context.contains("N"))
                expectedContext = context;
        }
        return expectedContext;
    }

    private static String stringFrom(byte[] array) {
        String s = "";
        for (byte value : array)
            s += (char) value;
        return s;
    }

}
