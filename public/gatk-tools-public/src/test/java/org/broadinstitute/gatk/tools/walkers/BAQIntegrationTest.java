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

package org.broadinstitute.gatk.tools.walkers;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class BAQIntegrationTest extends WalkerTest {
    private final static String baseCommand = "-T PrintReads -R " + b36KGReference +
            " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
            " -o %s" +
            " -L 1:10,000,000-10,100,000";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing BAQ with SLX, 454, and SOLID data
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testPrintReadsNoBAQ() {
        WalkerTestSpec spec = new WalkerTestSpec( baseCommand +" -baq OFF",  1, Arrays.asList("e33187ca383c7f5c75c5d547ec79e1cb"));
        executeTest(String.format("testPrintReadsNoBAQ"), spec);
    }

    @Test
    public void testPrintReadsRecalBAQ() {
        WalkerTestSpec spec = new WalkerTestSpec( baseCommand +" -baq RECALCULATE",  1, Arrays.asList("a25043edfbfa4f21a13cc21064b460df"));
        executeTest(String.format("testPrintReadsRecalBAQ"), spec);
    }
}
