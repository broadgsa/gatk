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
        return "-T VariantSelect -o %s -R " + b36KGReference;
    }


    @Test
    public void testVCFSelect1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("8fd97a99174483920f84aaf54eeb9dd9"));
        executeTest("testVCFSelect1", spec);
    }

    @Test
    public void testVCFSelect2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -L 1:10001290-10048590 ", 1,
                Arrays.asList("ecb61d798ff5e0b003d2986500bf4462"));
        executeTest("testVCFSelect2", spec);
    }

    @Test
    public void testVCFSelectOr() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("7c23412a0abe28057a78cd532b752019"));
        executeTest("testVCFSelectOr", spec);
    }

    @Test
    public void testVCFSelectAnd() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6 && AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("71cc18348a9d5de7270613bedb79880a"));
        executeTest("testVCFSelectAnd", spec);
    }
}