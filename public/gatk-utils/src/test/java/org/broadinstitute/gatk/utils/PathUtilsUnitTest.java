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
import org.testng.Assert;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


import java.io.File;

public class PathUtilsUnitTest extends BaseTest {
    @BeforeClass
    public void init() { }

    /**
     * Tests that we can successfully refresh a volume
     */
    @Test
    public void testRefreshVolume() {
        logger.warn("Executing testRefreshVolume");

        Assert.assertTrue(successfullyRefreshedVolume(System.getProperty("java.io.tmpdir")));
        Assert.assertFalse(successfullyRefreshedVolume("/a/made/up/file.txt"));
    }

    private boolean successfullyRefreshedVolume(String filename) {
        boolean result = true;

        try {
            PathUtils.refreshVolume(new File(filename));
        } catch (ReviewedGATKException e) {
            result = false;
        }

        logger.warn(filename + " is accessible : " + result);

        return result;
    }
}
