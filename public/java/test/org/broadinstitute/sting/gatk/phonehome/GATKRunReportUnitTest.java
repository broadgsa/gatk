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

package org.broadinstitute.sting.gatk.phonehome;

import junit.framework.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.testng.annotations.Test;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class GATKRunReportUnitTest extends BaseTest {
    @Test
    public void testAccessKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSAccessKey(), "c0f0afa1ff5ba41d9bf216cfcdbf26bf");
    }

    @Test
    public void testSecretKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSSecretKey(), "db2f13b3a7c98ad24e28783733ec4a62");
    }

    private void testAWSKey(final String accessKey, final String expectedMD5) throws Exception {
        Assert.assertNotNull(accessKey, "AccessKey should not be null");
        final String actualmd5 = Utils.calcMD5(accessKey);
        Assert.assertEquals(actualmd5, expectedMD5);
    }
}
