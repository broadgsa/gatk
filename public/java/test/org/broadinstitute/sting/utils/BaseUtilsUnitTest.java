package org.broadinstitute.sting.utils;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeClass;


public class BaseUtilsUnitTest extends BaseTest {
    @BeforeClass
    public void init() { }

    @Test
    public void testMostFrequentBaseFraction() {
        logger.warn("Executing testMostFrequentBaseFraction");
        
        compareFrequentBaseFractionToExpected("AAAAA", 1.0);
        compareFrequentBaseFractionToExpected("ACCG", 0.5);
        compareFrequentBaseFractionToExpected("ACCCCTTTTG", 4.0/10.0);
    }

    private void compareFrequentBaseFractionToExpected(String sequence, double expected) {
        double fraction = BaseUtils.mostFrequentBaseFraction(sequence.getBytes());
        Assert.assertTrue(MathUtils.compareDoubles(fraction, expected) == 0);
    }

    @Test
    public void testTransitionTransversion() {
        logger.warn("Executing testTransitionTransversion");

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
    }

    @Test
    public void testReverseComplementString() {
        logger.warn("Executing testReverseComplementString");

        compareRCStringToExpected("ACGGT", "ACCGT");
        compareRCStringToExpected("TCGTATATCTCGCTATATATATATAGCTCTAGTATA", "TATACTAGAGCTATATATATATAGCGAGATATACGA");
        compareRCStringToExpected("AAAN", "NTTT");
    }

    private void compareRCStringToExpected(String fw, String rcExp) {
        String rcObs = BaseUtils.simpleReverseComplement(fw);

        Assert.assertTrue(rcObs.equals(rcExp));
    }
}
