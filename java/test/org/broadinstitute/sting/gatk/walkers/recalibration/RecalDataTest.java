// our package
package org.broadinstitute.sting.gatk.walkers.recalibration;


// the imports for unit testing.

import org.junit.Assert;
import org.junit.Test;
import org.junit.Before;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.QualityUtils;
import java.util.HashSet;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Basic unit test for RecalData
 */
public class RecalDataTest extends BaseTest {
    RecalData datum, starDatum;

    /**
     * Tests that we got a string parameter in correctly
     */
    @Before
    public void before() {
        datum = new RecalData(2, 3, "0", "AA");
        datum.B = datum.N = 1;
        starDatum = new RecalData(2, 3, "0", "**");
        starDatum.B = starDatum.N = 1;
    }

    @Test
    public void testBasic() {
        logger.warn("Executing testIsBetween");

        Assert.assertTrue(datum.N == 1);
        Assert.assertTrue(datum.B == 1);
        Assert.assertTrue(datum.pos == 2);
        Assert.assertTrue(datum.qual == 3);
        Assert.assertTrue(datum.readGroup.equals("0"));
        Assert.assertTrue(datum.dinuc.equals("AA"));

        Assert.assertTrue(starDatum.N == 1);
        Assert.assertTrue(starDatum.B == 1);
        Assert.assertTrue(starDatum.pos == 2);
        Assert.assertTrue(starDatum.qual == 3);
        Assert.assertTrue(starDatum.readGroup.equals("0"));
        Assert.assertTrue(starDatum.dinuc.equals("**"));
    }

    @Test
    public void testInc() {
        logger.warn("Executing testInc");

        datum.inc(1L, 0L);
        Assert.assertTrue(datum.N == 2);
        Assert.assertTrue(datum.B == 1);

        datum.inc('A', 'A');
        Assert.assertTrue(datum.N == 3);
        Assert.assertTrue(datum.B == 1);

        datum.inc('A', 'C');
        Assert.assertTrue(datum.N == 4);
        Assert.assertTrue(datum.B == 2);
    }


    @Test
    public void testEmpQual() {
        logger.warn("Executing testEmpQual");

        datum.B = 1;
        datum.N = 1;
        Assert.assertEquals(datum.empiricalQualDouble(), 0.0, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  1);

        datum.B = 1;
        datum.N = 2;
        Assert.assertEquals(datum.empiricalQualDouble(), 3.0103, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  3);

        datum.B = 1;
        datum.N = 3;
        Assert.assertEquals(datum.empiricalQualDouble(),  4.771213, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  5);

        datum.B = 1;
        datum.B = 2;
        Assert.assertEquals(datum.empiricalQualDouble(),  1.760913, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  2);

        datum.B = 1;
        datum.N = 10;
        Assert.assertEquals(datum.empiricalQualDouble(),  10.0, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  10);

        datum.B = 1;
        datum.N = 100;
        Assert.assertEquals(datum.empiricalQualDouble(),  20.0, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  20);

        datum.B = 1;
        datum.N = 1000;
        Assert.assertEquals(datum.empiricalQualDouble(),  30.0, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  30);

        datum.B = 0;
        datum.N = 1000;
        Assert.assertEquals(datum.empiricalQualDouble(),  QualityUtils.MAX_REASONABLE_Q_SCORE, 1e-5);
        Assert.assertEquals(datum.empiricalQualByte(),  QualityUtils.MAX_REASONABLE_Q_SCORE);
    }


    public void testtoCSVString() {
        logger.warn("Executing testtoCSVString");

        Assert.assertEquals(datum.toCSVString(false), "0,AA,3,2,1,1,0");
        Assert.assertEquals(datum.toCSVString(true), "0,AA,3,*,1,1,0");
    }


    public void testFromCSVString() {
        logger.warn("Executing testFromCSVString");

        Assert.assertEquals(RecalData.fromCSVString("0,AA,3,2,1,1,0").toCSVString(false), datum.toCSVString(false));
        Assert.assertEquals(RecalData.fromCSVString("0,AA,3,*,1,1,0").toCSVString(false), datum.toCSVString(true));
        Assert.assertEquals(RecalData.fromCSVString("0,**,3,2,1,1,0").toCSVString(false), starDatum.toCSVString(false));
        Assert.assertEquals(RecalData.fromCSVString("0,**,3,*,1,1,0").toCSVString(false), starDatum.toCSVString(true));
    }

    public void testDinucIndex() {
        logger.warn("Executing testDinucIndex");

        HashSet<Integer> indices = new HashSet<Integer>();
        HashSet<Byte> unknownBytes = new HashSet<Byte>();
        byte bases[] = {'A', 'C', 'G', 'T', '*', 'N'};
        unknownBytes.add((byte)'*');
        unknownBytes.add((byte)'N');

        for ( int i = 0; i < bases.length; i++ ) {
            for ( int j = 0; j < bases.length; j++ ) {
                byte[] bp = {bases[i], bases[j]};
                String s = new String(bp);
                int index = RecalData.dinucIndex(s);
                indices.add(index);
                Assert.assertEquals(index, RecalData.dinucIndex(bases[i], bases[j]));
                if ( index != -1 ) {
                    Assert.assertEquals(RecalData.dinucIndex2bases(index), s);
                } else {
                    Assert.assertTrue(unknownBytes.contains(bp[0]) || unknownBytes.contains(bp[1]) );
                }
            }
        }
        Assert.assertEquals(indices.size(), 17);
    }
}