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

package org.broadinstitute.gatk.utils;

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
 * UnitTests for the AutoFormatting
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class AutoFormattingTimeUnitTest extends BaseTest {
    @DataProvider(name = "AutoFormattingTimeUnitSelection")
    public Object[][] makeTimeData() {
        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{TimeUnit.SECONDS.toNanos(10), "s"});
        tests.add(new Object[]{TimeUnit.MINUTES.toNanos(10), "m"});
        tests.add(new Object[]{TimeUnit.HOURS.toNanos(10), "h"});
        tests.add(new Object[]{TimeUnit.DAYS.toNanos(10), "d"});
        tests.add(new Object[]{TimeUnit.DAYS.toNanos(1000), "w"});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AutoFormattingTimeUnitSelection")
    public void testUnitSelection(final long nano, final String expectedUnit) throws InterruptedException {
        final AutoFormattingTime time = new AutoFormattingTime(nano);
        testBasic(time, nano, time.getWidth(), time.getPrecision());
        Assert.assertTrue(time.toString().endsWith(expectedUnit), "TimeUnit " + time.toString() + " didn't contain expected time unit " + expectedUnit);
    }

    @Test(dataProvider = "AutoFormattingTimeUnitSelection")
    public void testSecondsAsDouble(final long nano, final String expectedUnit) throws InterruptedException {
        final double inSeconds = nano * 1e-9;
        final long nanoFromSeconds = (long)(inSeconds * 1e9);
        final AutoFormattingTime time = new AutoFormattingTime(inSeconds);
        testBasic(time, nanoFromSeconds, time.getWidth(), time.getPrecision());
    }

    @DataProvider(name = "AutoFormattingTimeWidthAndPrecision")
    public Object[][] makeTimeWidthAndPrecision() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final int width : Arrays.asList(-1, 1, 2, 6, 20) ) {
            for ( final int precision : Arrays.asList(1, 2) ) {
                tests.add(new Object[]{100.123456 * 1e9, width, precision});
                tests.add(new Object[]{0.123456 * 1e9, width, precision});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AutoFormattingTimeWidthAndPrecision")
    public void testWidthAndPrecision(final double inSeconds, final int width, final int precision) throws InterruptedException {
        final AutoFormattingTime time = new AutoFormattingTime(inSeconds, width, precision);
        final long nanoFromSeconds = (long)(inSeconds * 1e9);
        testBasic(time, nanoFromSeconds, width, precision);
        final Matcher match = matchToString(time);
        match.matches();
        final String widthString = match.group(1);
        final String precisionString = match.group(2);
        if ( width != -1 ) {
            final int actualWidth = widthString.length() + 1 + precisionString.length();
            Assert.assertTrue(actualWidth >= width, "width string '" + widthString + "' not >= the expected width " + width);
        }
        Assert.assertEquals(precisionString.length(), precision, "precision string '" + precisionString + "' not the expected precision " + precision);
    }

    private static Matcher matchToString(final AutoFormattingTime time) {
        Pattern pattern = Pattern.compile("(\\s*\\d*)\\.(\\d*) \\w");
        return pattern.matcher(time.toString());
    }

    private static void testBasic(final AutoFormattingTime aft, final long nano, final int expectedWidth, final int expectedPrecision) {
        Assert.assertEquals(aft.getTimeInNanoSeconds(), nano);
        assertEqualsDoubleSmart(aft.getTimeInSeconds(), nano * 1e-9, 1e-3, "Time in seconds not within tolerance of nanoSeconds");
        Assert.assertEquals(aft.getWidth(), expectedWidth);
        Assert.assertEquals(aft.getPrecision(), expectedPrecision);
        Assert.assertNotNull(aft.toString(), "TimeUnit toString returned null");
        final Matcher match = matchToString(aft);
        Assert.assertTrue(match.matches(), "toString " + aft.toString() + " doesn't match our expected format");
    }
}
