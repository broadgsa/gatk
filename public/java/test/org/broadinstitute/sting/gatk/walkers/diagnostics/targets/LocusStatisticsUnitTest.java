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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Set;

public class LocusStatisticsUnitTest /*extends BaseTest*/ {

    @Test(dataProvider = "StatusTestValues")
    public void testCallableStatuses(int coverage, int rawCoverage, CallableStatus status) {
        // The min Coverage threshold is 10, the max is 100
        ThresHolder thresholds = new ThresHolder(20, 20, 10, 100, 20, 50, 0.5, 0.2, 0.5, 0.2, 0.2, 0.5);
        Set<CallableStatus> statuses = new LocusStatistics(coverage, rawCoverage).callableStatuses(thresholds);
        // Check to make sure the status provides matches the actual
        Assert.assertTrue((status == null) ? statuses.isEmpty() : (statuses.contains(status) && statuses.size() == 1));

    }

    @DataProvider(name = "StatusTestValues")
    public Object[][] getStatusTestValues() {
        return new Object[][]{
                new Object[]{100, 100, null},
                new Object[]{100, 101, null},
                new Object[]{101, 101, CallableStatus.EXCESSIVE_COVERAGE},
                new Object[]{10, 101, null},
                new Object[]{9, 101, CallableStatus.POOR_QUALITY},
                new Object[]{9, 10, CallableStatus.POOR_QUALITY},
                new Object[]{9, 9, CallableStatus.LOW_COVERAGE},
                new Object[]{0, 0, CallableStatus.COVERAGE_GAPS},
                new Object[]{0, 9, CallableStatus.LOW_COVERAGE},
                new Object[]{0, 101, CallableStatus.POOR_QUALITY},
                new Object[]{10, Integer.MAX_VALUE, null},
                new Object[]{Integer.MAX_VALUE, Integer.MAX_VALUE, CallableStatus.EXCESSIVE_COVERAGE},
        };
    }

}
