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

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class CallableLociIntegrationTest extends WalkerTest {
    final static String commonArgs     = "-R " + b36KGReference + " -T CallableLoci -I " + validationDataLocation + "/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s";

    final static String SUMMARY_MD5 = "a6f5963669f19d9d137ced87d65834b0";

    @Test
    public void testCallableLociWalkerBed() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("9b4ffea1dbcfefadeb1c9fa74b0e0e59", SUMMARY_MD5));
        executeTest("formatBed", spec);
    }

    @Test
    public void testCallableLociWalkerPerBase() {
        String gatk_args = commonArgs + " -format STATE_PER_BASE -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("d6505e489899e80c08a7168777f6e07b", SUMMARY_MD5));
        executeTest("format_state_per_base", spec);
    }
    
    @Test
    public void testCallableLociWalker2() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-10,000,100 -L 1:10,000,110-10,000,120 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("330f476085533db92a9dbdb3a127c041", "d287510eac04acf5a56f5cde2cba0e4a"));
        executeTest("formatBed by interval", spec);
    }

    @Test
    public void testCallableLociWalker3() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-11,000,000 -minDepth 10 -maxDepth 100 --minBaseQuality 10 --minMappingQuality 20 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("7f79ad8195c4161060463eeb21d2bb11", "7ee269e5f4581a924529a356cc806e55"));
        executeTest("formatBed lots of arguments", spec);
    }
}
