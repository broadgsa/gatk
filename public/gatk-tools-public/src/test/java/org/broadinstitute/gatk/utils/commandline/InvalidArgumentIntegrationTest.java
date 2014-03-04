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

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 8/31/12
 * Time: 11:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class InvalidArgumentIntegrationTest extends WalkerTest {
    private static final String callsB36  = BaseTest.validationDataLocation + "lowpass.N3.chr1.raw.vcf";

    private WalkerTest.WalkerTestSpec baseTest(String flag, String arg, Class exeption) {
        return new WalkerTest.WalkerTestSpec("-T VariantsToTable -M 10 --variant:vcf "
                + callsB36 + " -F POS,CHROM -R "
                + b36KGReference +  " -o %s " + flag + " " + arg,
                1, exeption);

    }

    @Test
    public void testUnknownReadFilter() {
        executeTest("UnknownReadFilter",baseTest("-rf","TestUnknownReadFilter", UserException.MalformedReadFilterException.class));
    }

    @Test
    public void testMalformedWalkerArgs() {
        executeTest("MalformedWalkerArgs",
                new WalkerTest.WalkerTestSpec("-T UnknownWalkerName -M 10 --variant:vcf "
                + callsB36 + " -F POS,CHROM -R "
                + b36KGReference +  " -o %s ",
                1, UserException.MalformedWalkerArgumentsException.class));
    }
}
