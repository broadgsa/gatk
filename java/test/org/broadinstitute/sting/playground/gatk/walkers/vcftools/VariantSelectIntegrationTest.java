/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantSelectIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantSelect -o %s -R " + oneKGLocation + "reference/human_b36_both.fasta";
    }


    @Test
    public void testVCFSelect1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("8b358e0cfa35de022a37360a6f28a839"));
        executeTest("testVCFSelect1", spec);
    }

    @Test
    public void testVCFSelect2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -L 1:10001290-10048590 ", 1,
                Arrays.asList("8e991b9d6d610c8f89c8557756fc8e34"));
        executeTest("testVCFSelect2", spec);
    }

    @Test
    public void testVCFSelectOr() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("7bd064c8d8bf5389fcd0b78a7c32b599"));
        executeTest("testVCFSelectOr", spec);
    }

    @Test
    public void testVCFSelectAnd() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6 && AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("5af565836fa926feaa130715b93188bc"));
        executeTest("testVCFSelectAnd", spec);
    }
}