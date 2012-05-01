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

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/21/12
 */

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Basic unit test for QualityUtils class
 */
public class QualityUtilsUnitTest extends BaseTest {
    @BeforeClass
    public void init() {
    }

    @Test
    public void testQualCaches() {
        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 20), 0.01, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 20), -2.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 20), 0.99, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 20), -0.0043648054, 1e-6);

        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 30), 0.001, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 30), -3.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 30), 0.999, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 30), -0.000434511774, 1e-6);

        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 40), 0.0001, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 40), -4.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 40), 0.9999, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 40), -4.34316198e-5, 1e-6);
    }
}