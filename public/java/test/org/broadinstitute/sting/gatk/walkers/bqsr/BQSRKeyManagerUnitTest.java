package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Mauricio Carneiro
 * @since 3/7/12
 */
public class BQSRKeyManagerUnitTest {
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
    }

    @Test(enabled = false)
    public void testCombineBitSets() {
        final int nRequired = 2;
        final ArrayList<Covariate> covariates = new ArrayList<Covariate>();
        covariates.add(new ReadGroupCovariate());
        covariates.add(new QualityScoreCovariate());
        covariates.add(new CycleCovariate());
        covariates.add(new ContextCovariate());
        createReadAndTest(covariates, nRequired);
    }
    
    @Test(enabled = true)
    public void testOnlyRequiredCovariates() {
        final int nRequired = 2;
        final ArrayList<Covariate> covariates = new ArrayList<Covariate>(2);
        covariates.add(new ReadGroupCovariate());
        covariates.add(new QualityScoreCovariate());
        createReadAndTest(covariates, nRequired);
    }

    @Test(enabled = true)
    public void testOnlyOneCovariate() {
        final int nRequired = 1;
        final ArrayList<Covariate> covariates = new ArrayList<Covariate>(2);
        covariates.add(new ReadGroupCovariate());
        createReadAndTest(covariates, nRequired);
    }

    @Test(enabled = false)
    public void testOneCovariateWithOptionalCovariates() {
        final int nRequired = 1;
        final ArrayList<Covariate> covariates = new ArrayList<Covariate>(4);
        covariates.add(new ReadGroupCovariate());
        covariates.add(new QualityScoreCovariate());
        covariates.add(new CycleCovariate());
        covariates.add(new ContextCovariate());
        createReadAndTest(covariates, nRequired);
    }

    private void createReadAndTest(List<Covariate> covariates, int nRequired) {
        int readLength = 1000;
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(ReadUtils.createRandomReadBases(readLength, true), ReadUtils.createRandomReadQuals(readLength), readLength + "M");
        read.setReadGroup(new GATKSAMReadGroupRecord("ID"));
        read.getReadGroup().setPlatform("illumina");

        runTestOnRead(read, covariates, nRequired);
        read.setReadNegativeStrandFlag(true);
        runTestOnRead(read, covariates, nRequired);
        read.setReadPairedFlag(true);
        read.setSecondOfPairFlag(true);
        runTestOnRead(read, covariates, nRequired);
        read.setReadNegativeStrandFlag(false);
        runTestOnRead(read, covariates, nRequired);
    }

    private void runTestOnRead(GATKSAMRecord read, List<Covariate> covariateList, int nRequired) {
        final long[][][] covariateKeys = new long[covariateList.size()][EventType.values().length][];
        int i = 0;
        for (Covariate cov : covariateList) {
            cov.initialize(RAC);
            CovariateValues covValues = cov.getValues(read);
            covariateKeys[i][EventType.BASE_SUBSTITUTION.index] = covValues.getMismatches();
            covariateKeys[i][EventType.BASE_INSERTION.index] = covValues.getInsertions();
            covariateKeys[i][EventType.BASE_DELETION.index] = covValues.getDeletions();
            i++;
        }
        List<Covariate> requiredCovariates = new LinkedList<Covariate>();
        List<Covariate> optionalCovariates = new LinkedList<Covariate>();
        
        for (int j=0; j<nRequired; j++)
            requiredCovariates.add(covariateList.get(j));
        for (int j=nRequired; j<covariateList.size(); j++)
            optionalCovariates.add(covariateList.get(j));
            
        BQSRKeyManager keyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);

        for (int l = 0; l < read.getReadLength(); l++) {
            for (EventType eventType : EventType.values()) {
                long[] keySet = new long[covariateList.size()];
                Object[] expectedRequired = new Object[covariateList.size()];
                Object[] expectedCovariate = new Object[covariateList.size()];

                for (int j = 0; j < covariateList.size(); j++) {
                    keySet[j] = covariateKeys[j][eventType.index][l];

                    if (j < nRequired)
                        expectedRequired[j] = covariateList.get(j).formatKey(keySet[j]);
                    else
                        expectedCovariate[j - nRequired] = covariateList.get(j).formatKey(keySet[j]);
                }

                if (optionalCovariates.size() == 0) {
                    final long masterKey = keyManager.createMasterKey(keySet, eventType, -1);
                    testKeys(keyManager, masterKey, nRequired, optionalCovariates, expectedRequired, expectedCovariate, eventType, -1);
                } else {
                    for (int j = 0; j < optionalCovariates.size(); j++) {
                        final long masterKey = keyManager.createMasterKey(keySet, eventType, j);
                        testKeys(keyManager, masterKey, nRequired, optionalCovariates, expectedRequired, expectedCovariate, eventType, j);
                    }
                }
            }
        }
    }

    private void testKeys(final BQSRKeyManager keyManager, final long key, final int nRequired, final List<Covariate> optionalCovariates,
                          final Object[] expectedRequired, final Object[] expectedCovariate, final EventType eventType, final int index) {

        Object[] actual = keyManager.keySetFrom(key).toArray();

        // Build the expected array
        Object[] expected = new Object[nRequired + (optionalCovariates.size() > 0 ? 3 : 1)];
        System.arraycopy(expectedRequired, 0, expected, 0, nRequired);
        if (optionalCovariates.size() > 0) {
            expected[expected.length-3] = expectedCovariate[index];
            expected[expected.length-2] = optionalCovariates.get(index).getClass().getSimpleName().split("Covariate")[0];
        }
        expected[expected.length-1] = eventType;

//                    System.out.println("Actual  : " + Utils.join(",", Arrays.asList(actual)));
//                    System.out.println("Expected: " + Utils.join(",", Arrays.asList(expected)));
//                    System.out.println();

        for (int k = 0; k < expected.length; k++)
            Assert.assertEquals(actual[k], expected[k]);
    }
}
