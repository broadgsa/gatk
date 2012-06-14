package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.clipping.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
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

    @Test(enabled = true)
    public void testSimpleContexts() {
        GATKSAMRecord read = ReadUtils.createRandomRead(1000);
        GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        CovariateValues values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(values.getInsertions(), RAC.INSERTIONS_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(values.getDeletions(),  RAC.DELETIONS_CONTEXT_SIZE,  clippedRead, covariate);
    }

    public static void verifyCovariateArray(long[] values, int contextSize, GATKSAMRecord read, Covariate contextCovariate) {
        for (int i = 0; i < values.length; i++)
            Assert.assertEquals(contextCovariate.formatKey(values[i]), expectedContext(read, i, contextSize));

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
