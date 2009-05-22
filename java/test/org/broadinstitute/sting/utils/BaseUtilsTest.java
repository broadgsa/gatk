package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.BaseTest;
import org.junit.Test;
import org.junit.BeforeClass;
import org.junit.Assert;

public class BaseUtilsTest extends BaseTest {
    @BeforeClass
    public static void init() { }

    @Test
    public void testMostFrequentBaseFraction() {
        logger.warn("Executing testMostFrequentBaseFraction");
        
        compareFrequentBaseFractionToExpected("AAAAA", 1.0);
        compareFrequentBaseFractionToExpected("ACCG", 0.5);
        compareFrequentBaseFractionToExpected("ACCCCTTTTG", 4.0/10.0);
    }

    @Test
    public void testTransitionTransversion() {
        logger.warn("Executing testTransitionTransversion");

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'G' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'C', 'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'C', 'T' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'C', 'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'T', 'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'T', 'C' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'T', 'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'G', 'A' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'G', 'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'G', 'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'a', 'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'a', 'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 't' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'A', 'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'a', 't' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( 'a', 'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
    }

    private void compareFrequentBaseFractionToExpected(String sequence, double expected) {
        double fraction = BaseUtils.mostFrequentBaseFraction(sequence.getBytes());
        Assert.assertTrue(MathUtils.compareDoubles(fraction, expected) == 0);
    }
}
