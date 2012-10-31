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


public class RecalDatumUnitTest extends BaseTest {

    // --------------------------------------------------------------------------------
    //
    // merge case Provider
    //
    // --------------------------------------------------------------------------------

    private class RecalDatumTestProvider extends TestDataProvider {
        int exError, exTotal, reportedQual;

        private RecalDatumTestProvider(int E, int N, int reportedQual) {
            super(RecalDatumTestProvider.class);

            this.exError = E;
            this.exTotal = N;
            this.reportedQual = reportedQual;
        }

        public double getErrorRate() {
            return (exError + 1) / (1.0 * (exTotal + 2));
        }

        public double getErrorRatePhredScaled() {
            return QualityUtils.phredScaleErrorRate(getErrorRate());
        }

        public int getReportedQual() {
            return reportedQual;
        }

        public RecalDatum makeRecalDatum() {
            return new RecalDatum(exTotal, exError, (byte)getReportedQual());
        }

        @Override
        public String toString() {
            return String.format("exError=%d, exTotal=%d, reportedQual=%d", exError, exTotal, reportedQual);
        }
    }

    @DataProvider(name = "RecalDatumTestProvider")
    public Object[][] makeRecalDatumTestProvider() {
        for ( int E : Arrays.asList(1, 10, 100, 1000, 10000) )
            for ( int N : Arrays.asList(10, 100, 1000, 10000, 100000, 1000000) )
                for ( int reportedQual : Arrays.asList(10, 20) )
                    if ( E <= N )
                        new RecalDatumTestProvider(E, N, reportedQual);
        return RecalDatumTestProvider.getTests(RecalDatumTestProvider.class);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumBasics(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }

    private static void assertBasicFeaturesOfRecalDatum(final RecalDatum datum, final RecalDatumTestProvider cfg) {
        Assert.assertEquals(datum.getNumMismatches(), cfg.exError, 1E-6);
        Assert.assertEquals(datum.getNumObservations(), cfg.exTotal, 1E-6);
        if ( cfg.getReportedQual() != -1 )
            Assert.assertEquals(datum.getEstimatedQReportedAsByte(), cfg.getReportedQual());
        BaseTest.assertEqualsDoubleSmart(datum.getEmpiricalQuality(), cfg.getErrorRatePhredScaled());
        BaseTest.assertEqualsDoubleSmart(datum.getEmpiricalErrorRate(), cfg.getErrorRate());
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumCopyAndCombine(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        final RecalDatum copy = new RecalDatum(datum);
        assertBasicFeaturesOfRecalDatum(copy, cfg);

        RecalDatumTestProvider combinedCfg = new RecalDatumTestProvider(cfg.exError * 2, cfg.exTotal * 2, cfg.reportedQual);
        copy.combine(datum);
        assertBasicFeaturesOfRecalDatum(copy, combinedCfg);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumModification(RecalDatumTestProvider cfg) {
        RecalDatum datum = cfg.makeRecalDatum();
        datum.setEmpiricalQuality(10.1);
        Assert.assertEquals(datum.getEmpiricalQuality(), 10.1);

        datum.setEstimatedQReported(10.1);
        Assert.assertEquals(datum.getEstimatedQReported(), 10.1);
        Assert.assertEquals(datum.getEstimatedQReportedAsByte(), 10);

        datum = cfg.makeRecalDatum();
        cfg.exTotal = 100000;
        datum.setNumObservations(cfg.exTotal);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        cfg.exError = 1000;
        datum.setNumMismatches(cfg.exError);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(true);
        cfg.exError++;
        cfg.exTotal++;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(10, 5);
        cfg.exError += 5;
        cfg.exTotal += 10;
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }
}