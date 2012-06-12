package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * @author carneiro
 * @since 4/21/12
 */
public class ReadCovariatesUnitTest {

    @Test(enabled = false)
    public void testCovariateGeneration() {
        final String RGID = "id";
        final int length = 10;
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        GATKSAMRecord read = ReadUtils.createRandomRead(length, false);
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(RGID);
        rg.setPlatform("illumina");
        read.setReadGroup(rg);
        final byte[] mQuals = read.getBaseQualities(EventType.BASE_SUBSTITUTION);
        final byte[] iQuals = read.getBaseQualities(EventType.BASE_INSERTION);
        final byte[] dQuals = read.getBaseQualities(EventType.BASE_DELETION);

        ReadGroupCovariate rgCov = new ReadGroupCovariate();
        QualityScoreCovariate qsCov = new QualityScoreCovariate();
        ContextCovariate coCov = new ContextCovariate();
        CycleCovariate cyCov = new CycleCovariate();

        rgCov.initialize(RAC);
        qsCov.initialize(RAC);
        coCov.initialize(RAC);
        cyCov.initialize(RAC);

        List<Covariate> requestedCovariates = new ArrayList<Covariate>(4);
        requestedCovariates.add(rgCov);
        requestedCovariates.add(qsCov);
        requestedCovariates.add(coCov);
        requestedCovariates.add(cyCov);

        ReadCovariates rc = RecalDataManager.computeCovariates(read, requestedCovariates);

        // check that the length is correct
        Assert.assertEquals(rc.getMismatchesKeySet().length, length);
        Assert.assertEquals(rc.getInsertionsKeySet().length, length);
        Assert.assertEquals(rc.getDeletionsKeySet().length, length);

        for (int i = 0; i < length; i++) {
            // check that read group is always the same
            Assert.assertEquals(rgCov.formatKey(rc.getMismatchesKeySet(i)[0]), RGID);
            Assert.assertEquals(rgCov.formatKey(rc.getInsertionsKeySet(i)[0]), RGID);
            Assert.assertEquals(rgCov.formatKey(rc.getDeletionsKeySet(i)[0]),  RGID);

            // check quality score
            Assert.assertEquals(qsCov.formatKey(rc.getMismatchesKeySet(i)[1]), "" + mQuals[i]);
            Assert.assertEquals(qsCov.formatKey(rc.getInsertionsKeySet(i)[1]), "" + iQuals[i]);
            Assert.assertEquals(qsCov.formatKey(rc.getDeletionsKeySet(i)[1]),  "" + dQuals[i]);

            // check context
            Assert.assertEquals(coCov.formatKey(rc.getMismatchesKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.MISMATCHES_CONTEXT_SIZE));
            Assert.assertEquals(coCov.formatKey(rc.getInsertionsKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.INSERTIONS_CONTEXT_SIZE));
            Assert.assertEquals(coCov.formatKey(rc.getDeletionsKeySet(i)[2]),  ContextCovariateUnitTest.expectedContext(read, i, RAC.DELETIONS_CONTEXT_SIZE));

            // check cycle
            Assert.assertEquals(cyCov.formatKey(rc.getMismatchesKeySet(i)[3]), "" + (i+1));
            Assert.assertEquals(cyCov.formatKey(rc.getInsertionsKeySet(i)[3]), "" + (i+1));
            Assert.assertEquals(cyCov.formatKey(rc.getDeletionsKeySet(i)[3]),  "" + (i+1));
        }

    }

}
