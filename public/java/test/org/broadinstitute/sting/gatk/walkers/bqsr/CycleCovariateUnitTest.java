package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.BitSet;
import java.util.Random;

/**
 * @author Mauricio Carneiro
 * @since 3/1/12
 */
public class CycleCovariateUnitTest {
    CycleCovariate covariate;
    RecalibrationArgumentCollection RAC;
    Random random;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new CycleCovariate();
        random = GenomeAnalysisEngine.getRandomGenerator();
        covariate.initialize(RAC);
    }

    @Test(enabled = true)
    public void testSimpleCycles() {
        short readLength = 10;
        byte[] quals = ReadUtils.createRandomReadQuals(readLength);
        byte[] bbases = ReadUtils.createRandomReadBases(readLength, true);
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bbases, quals, bbases.length + "M");
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        CovariateValues values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), (short) 1, (short) 1);

        read.setReadNegativeStrandFlag(true);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), readLength, (short) -1);

        read.setReadPairedFlag(true);
        read.setSecondOfPairFlag(true);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), (short) -readLength, (short) 1);

        read.setReadNegativeStrandFlag(false);
        values = covariate.getValues(read);
        verifyCovariateArray(values.getMismatches(), (short) -1, (short) -1);

    }

    private void verifyCovariateArray(BitSet[] values, short init, short increment) {
        for (short i = 0; i < values.length; i++) {
            short actual = Short.decode(covariate.keyFromBitSet(values[i]));
            int expected = init + (increment * i);
            //            System.out.println(String.format("%d: %d, %d", i, actual, expected));
            Assert.assertEquals(actual, expected);
        }
    }

}
