/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// our package
package org.broadinstitute.sting.utils.recalibration;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class QualQuantizerUnitTest extends BaseTest {
    @BeforeSuite
    public void before() {

    }

    // --------------------------------------------------------------------------------
    //
    // merge case Provider
    //
    // --------------------------------------------------------------------------------

    private class QualIntervalTestProvider extends TestDataProvider {
        final QualQuantizer.QualInterval left, right;
        int exError, exTotal, exQual;
        double exErrorRate;

        private QualIntervalTestProvider(int leftE, int leftN, int rightE, int rightN, int exError, int exTotal) {
            super(QualIntervalTestProvider.class);

            QualQuantizer qq = new QualQuantizer(0);
            left = qq.new QualInterval(10, 10, leftN, leftE, 0);
            right = qq.new QualInterval(11, 11, rightN, rightE, 0);

            this.exError = exError;
            this.exTotal = exTotal;
            this.exErrorRate = (leftE + rightE + 1) / (1.0 * (leftN + rightN + 1));
            this.exQual = QualityUtils.probToQual(1-this.exErrorRate, 0);
        }
    }

    @DataProvider(name = "QualIntervalTestProvider")
    public Object[][] makeQualIntervalTestProvider() {
        new QualIntervalTestProvider(10, 100, 10, 1000, 20, 1100);
        new QualIntervalTestProvider(0, 100, 10, 900,   10, 1000);
        new QualIntervalTestProvider(10, 900, 0, 100,   10, 1000);
        new QualIntervalTestProvider(0, 0, 10, 100,     10, 100);
        new QualIntervalTestProvider(1, 10, 9, 90,      10, 100);
        new QualIntervalTestProvider(1, 10, 9, 100000,  10, 100010);
        new QualIntervalTestProvider(1, 10, 9, 1000000, 10,1000010);

        return QualIntervalTestProvider.getTests(QualIntervalTestProvider.class);
    }

    @Test(dataProvider = "QualIntervalTestProvider")
    public void testQualInterval(QualIntervalTestProvider cfg) {
        QualQuantizer.QualInterval merged = cfg.left.merge(cfg.right);
        Assert.assertEquals(merged.nErrors, cfg.exError);
        Assert.assertEquals(merged.nObservations, cfg.exTotal);
        Assert.assertEquals(merged.getErrorRate(), cfg.exErrorRate);
        Assert.assertEquals(merged.getQual(), cfg.exQual);
    }

    @Test
    public void testMinInterestingQual() {
        for ( int q = 0; q < 15; q++ ) {
            for ( int minQual = 0; minQual <= 10; minQual ++ ) {
                QualQuantizer qq = new QualQuantizer(minQual);
                QualQuantizer.QualInterval left = qq.new QualInterval(q, q, 100, 10, 0);
                QualQuantizer.QualInterval right = qq.new QualInterval(q+1, q+1, 1000, 100, 0);

                QualQuantizer.QualInterval merged = left.merge(right);
                boolean shouldBeFree = q+1 <= minQual;
                if ( shouldBeFree )
                    Assert.assertEquals(merged.getPenalty(), 0.0);
                else
                    Assert.assertTrue(merged.getPenalty() > 0.0);
            }
        }
    }


    // --------------------------------------------------------------------------------
    //
    // High-level case Provider
    //
    // --------------------------------------------------------------------------------

    private class QuantizerTestProvider extends TestDataProvider {
        final List<Long> nObservationsPerQual = new ArrayList<Long>();
        final int nLevels;
        final List<Integer> expectedMap;

        private QuantizerTestProvider(final List<Integer> nObservationsPerQual, final int nLevels, final List<Integer> expectedMap) {
            super(QuantizerTestProvider.class);

            for ( int x : nObservationsPerQual )
                this.nObservationsPerQual.add((long)x);
            this.nLevels = nLevels;
            this.expectedMap = expectedMap;
        }

        @Override
        public String toString() {
            return String.format("QQTest nLevels=%d nObs=[%s] map=[%s]",
                    nLevels, Utils.join(",", nObservationsPerQual), Utils.join(",", expectedMap));
        }
    }

    @DataProvider(name = "QuantizerTestProvider")
    public Object[][] makeQuantizerTestProvider() {
        List<Integer> allQ2 = Arrays.asList(0, 0, 1000, 0, 0);

        new QuantizerTestProvider(allQ2, 5, Arrays.asList(0, 1, 2, 3, 4));
        new QuantizerTestProvider(allQ2, 1, Arrays.asList(2, 2, 2, 2, 2));

        new QuantizerTestProvider(Arrays.asList(0, 0, 1000, 0,  1000), 2, Arrays.asList(2, 2, 2, 2, 4));
        new QuantizerTestProvider(Arrays.asList(0, 0, 1000, 1,  1000), 2, Arrays.asList(2, 2, 2, 4, 4));
        new QuantizerTestProvider(Arrays.asList(0, 0, 1000, 10, 1000), 2, Arrays.asList(2, 2, 2, 2, 4));

        return QuantizerTestProvider.getTests(QuantizerTestProvider.class);
    }

    @Test(dataProvider = "QuantizerTestProvider", enabled = true)
    public void testQuantizer(QuantizerTestProvider cfg) {
        QualQuantizer qq = new QualQuantizer(cfg.nObservationsPerQual, cfg.nLevels, 0);
        logger.warn("cfg: " + cfg);
        for ( int i = 0; i < cfg.expectedMap.size(); i++) {
            int expected = cfg.expectedMap.get(i);
            int observed = qq.originalToQuantizedMap.get(i);
            //logger.warn(String.format("  qq map: %s : %d => %d", i, expected, observed));
            Assert.assertEquals(observed, expected);
        }
    }
}