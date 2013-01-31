/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
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
    public void testConvertIUPACtoN() {

        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'A'}, false, false), new byte[]{'A', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'W', 'A', 'A'}, false, false), new byte[]{'N', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'M', 'A'}, false, false), new byte[]{'A', 'N', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'K'}, false, false), new byte[]{'A', 'A', 'N'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'M', 'M', 'M'}, false, false), new byte[]{'N', 'N', 'N'});
    }

    private void checkBytesAreEqual(final byte[] b1, final byte[] b2) {
        for ( int i = 0; i < b1.length; i++ )
            Assert.assertEquals(b1[i], b2[i]);
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
