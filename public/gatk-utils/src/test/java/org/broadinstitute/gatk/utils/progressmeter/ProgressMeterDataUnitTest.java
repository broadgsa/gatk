/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.progressmeter;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.AutoFormattingTime;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * UnitTests for the ProgressMeterData
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ProgressMeterDataUnitTest extends BaseTest {
    @Test
    public void testBasic() {
        Assert.assertEquals(new ProgressMeterData(1.0, 2, 3).getElapsedSeconds(), 1.0);
        Assert.assertEquals(new ProgressMeterData(1.0, 2, 3).getUnitsProcessed(), 2);
        Assert.assertEquals(new ProgressMeterData(1.0, 2, 3).getBpProcessed(), 3);
    }

    @Test
    public void testFraction() {
        final double TOL = 1e-3;
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 1).calculateFractionGenomeTargetCompleted(10), 0.1, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 2).calculateFractionGenomeTargetCompleted(10), 0.2, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 1).calculateFractionGenomeTargetCompleted(100), 0.01, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 2).calculateFractionGenomeTargetCompleted(100), 0.02, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 1).calculateFractionGenomeTargetCompleted(0), 1.0, TOL);
    }

    @Test
    public void testSecondsPerBP() {
        final double TOL = 1e-3;
        final long M = 1000000;
        Assert.assertEquals(new ProgressMeterData(1.0, 1, M).secondsPerMillionBP(), 1.0, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, M/10).secondsPerMillionBP(), 10.0, TOL);
        Assert.assertEquals(new ProgressMeterData(2.0, 1, M).secondsPerMillionBP(), 2.0, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 1, 0).secondsPerMillionBP(), 1e6, TOL);
    }

    @Test
    public void testSecondsPerElement() {
        final double TOL = 1e-3;
        final long M = 1000000;
        Assert.assertEquals(new ProgressMeterData(1.0, M, 1).secondsPerMillionElements(), 1.0, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, M/10, 1).secondsPerMillionElements(), 10.0, TOL);
        Assert.assertEquals(new ProgressMeterData(2.00, M, 1).secondsPerMillionElements(), 2.0, TOL);
        Assert.assertEquals(new ProgressMeterData(1.0, 0, 1).secondsPerMillionElements(), 1e6, TOL);
    }
}
