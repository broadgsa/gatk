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


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ProcessUtilsUnitTest extends BaseTest {
    @Test()
    public void testGoodCommand() {
        final String goodLs = "ls " + b37KGReference;
        final int result = ProcessUtils.runCommandAndWait(goodLs);
        Assert.assertEquals(result, 0, "ProcessUtils tells me that my command that should be good failed");
    }

    @Test()
    public void testGoodCommandWithBadArguments() {
        final String goodLs = "ls asdfhadsfhakdhsfakdhfalkdhfalkdhflakhdflakdhsf";
        final int result = ProcessUtils.runCommandAndWait(goodLs);
        Assert.assertFalse(result == 0, "ProcessUtils tells me that my command which had a bad argument and should have returned not zero, did in fact return properly");
    }

    @Test(expectedExceptions = ReviewedStingException.class)
    public void testBadCommand() {
        final String goodLs = "asfdadsfadsfa " + b37KGReference;
        final int result = ProcessUtils.runCommandAndWait(goodLs);
        Assert.fail("ProcessUtils should have excepted out but got result back of " + result);
    }
}