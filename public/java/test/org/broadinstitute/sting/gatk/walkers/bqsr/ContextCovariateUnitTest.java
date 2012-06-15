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
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);

        verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INSERTIONS_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getDeletionsKeySet(),  RAC.DELETIONS_CONTEXT_SIZE,  clippedRead, covariate);
    }

    public static void verifyCovariateArray(long[][] values, int contextSize, GATKSAMRecord read, Covariate contextCovariate) {
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
