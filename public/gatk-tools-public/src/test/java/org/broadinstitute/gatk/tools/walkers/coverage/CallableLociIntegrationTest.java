/*
* Copyright 2012-2016 Broad Institute, Inc.
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

    final static String SUMMARY_MD5 = "27aea9e76af554db5e886dacb0d3edd0";

    @Test
    public void testCallableLociWalkerBed() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("908a22c0bfe81e6f4c571de290e48d03", SUMMARY_MD5));
        executeTest("formatBed", spec);
    }

    @Test
    public void testCallableLociWalkerPerBase() {
        String gatk_args = commonArgs + " -format STATE_PER_BASE -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("11e7fdd39bb83eec44ca452faf0dc825", SUMMARY_MD5));
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
                Arrays.asList("5e5ee8850815810680f205a520d05a0f", "441b7c31bed5c30705f9565b708ecace"));
        executeTest("formatBed lots of arguments", spec);
    }
}
