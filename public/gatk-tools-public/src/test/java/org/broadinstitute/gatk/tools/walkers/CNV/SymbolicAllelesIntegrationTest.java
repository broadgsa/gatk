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

package org.broadinstitute.gatk.tools.walkers.CNV;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SymbolicAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF) {
        return "-T CombineVariants" +
                " -R " + reference +
                " --variant:vcf " + privateTestDir + VCF +
                " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
                " -genotypeMergeOptions REQUIRE_UNIQUE" +
                " -setKey null" +
                " -o %s" +
                " --no_cmdline_in_header";
    }

    @Test(enabled = true)
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_1.vcf"),
                1,
                Arrays.asList("a1de53ac340f4ca02367c40680628251"));
        executeTest("Test symbolic alleles", spec);
    }

    @Test(enabled = true)
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_2.vcf"),
                1,
                Arrays.asList("c8b294089832bb1a2c450b550318a471"));
        executeTest("Test symbolic alleles mixed in with non-symbolic alleles", spec);
    }
}
