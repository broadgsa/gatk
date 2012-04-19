package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.clipping.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.BitSet;

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
        verifyCovariateArray(values.getMismatches(), RAC.MISMATCHES_CONTEXT_SIZE, stringFrom(clippedRead.getReadBases()));
        verifyCovariateArray(values.getInsertions(), RAC.INSERTIONS_CONTEXT_SIZE, stringFrom(clippedRead.getReadBases()));
        verifyCovariateArray(values.getDeletions(),  RAC.DELETIONS_CONTEXT_SIZE,  stringFrom(clippedRead.getReadBases()));
    }

    private void verifyCovariateArray(BitSet[] values, int contextSize, String bases) {
        for (int i = 0; i < values.length; i++) {
            String expectedContext = null;
            if (i - contextSize + 1 >= 0) {
                String context = bases.substring(i - contextSize + 1, i + 1);
                if (!context.contains("N"))
                    expectedContext = context;
            }
            Assert.assertEquals(covariate.keyFromBitSet(values[i]), expectedContext);
        }
    }

    private String stringFrom(byte[] array) {
        String s = "";
        for (byte value : array)
            s += (char) value;
        return s;
    }

}
