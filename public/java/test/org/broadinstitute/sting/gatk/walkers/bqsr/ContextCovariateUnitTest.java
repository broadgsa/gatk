package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.BitSet;
import java.util.Random;

/**
 * Short one line description of the walker.
 *
 * <p>
 * [Long description of the walker]
 * </p>
 *
 *
 * <h2>Input</h2>
 * <p>
 * [Description of the Input]
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * [Description of the Output]
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T [walker name]
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 3/1/12
 */
public class ContextCovariateUnitTest {
    ContextCovariate covariate;
    RecalibrationArgumentCollection RAC;
    Random random; 

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ContextCovariate();
        random = GenomeAnalysisEngine.getRandomGenerator();
        covariate.initialize(RAC);

    }

    @Test(enabled = true)
    public void testSimpleContexts() {
        byte [] quals = createRandomReadQuals(101);
        byte [] bbases = createRandomReadBases(101);
        String bases = stringFrom(bbases);
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bbases, quals, bbases.length + "M");
        CovariateValues values = covariate.getValues(read);
        verifyCovariateArray((BitSet []) values.getMismatches(), RAC.MISMATCHES_CONTEXT_SIZE, bases);
        verifyCovariateArray((BitSet []) values.getInsertions(), RAC.INSERTIONS_CONTEXT_SIZE, bases);
        verifyCovariateArray((BitSet []) values.getDeletions(),  RAC.DELETIONS_CONTEXT_SIZE, bases);
    }
    
    private void verifyCovariateArray(BitSet[] values, int contextSize, String bases) {
        for (int i=0; i<values.length; i++) {
            if (i >= contextSize)
                Assert.assertEquals(MathUtils.dnaFrom(values[i]), bases.substring(i-contextSize, i));
            else
                Assert.assertNull(values[i]);
        }
    }

    private String stringFrom(byte [] array) {
        String s = "";
        for (byte value : array)
            s += (char) value;
        return s;
    }

    private byte [] createRandomReadQuals(int length) {
        byte [] quals = new byte[length];
        for (int i=0; i<length; i++)
            quals[i] = (byte) random.nextInt(50);
        return quals;
    }

    private byte [] createRandomReadBases(int length) {
        byte [] bases = new byte[length];
        for (int i=0; i<length; i++) {
            switch(random.nextInt(4)) {
                case 0: bases[i] = 'A'; break;
                case 1: bases[i] = 'C'; break;
                case 2: bases[i] = 'G'; break;
                case 3: bases[i] = 'T'; break;
            }
        }
        return bases;
    }
}
