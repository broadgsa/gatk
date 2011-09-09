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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ValidateVariantsIntegrationTest extends WalkerTest {

    public static String baseTestString(String file, String type) {
        return "-T ValidateVariants -R " + b36KGReference + " -L 1:10001292-10001303 --variant:vcf " + validationDataLocation + file + " --validationType " + type;
    }

    @Test
    public void testGoodFile() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleGood.vcf", "ALL"),
                0,
                Arrays.asList("d41d8cd98f00b204e9800998ecf8427e")
        );

        executeTest("test good file", spec);
    }
    
    @Test
    public void testBadRefBase1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "REF"),
                0,
                UserException.MalformedFile.class
        );

        executeTest("test bad ref base #1", spec);
    }

    @Test
    public void testBadRefBase2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad2.vcf", "REF"),
                0,
                UserException.MalformedFile.class
        );

        executeTest("test bad ref base #2", spec);
    }

    @Test
    public void testBadChrCount1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "CHR_COUNTS"),
                0,
                UserException.MalformedFile.class
        );

        executeTest("test bad chr counts #1", spec);
    }

    @Test
    public void testBadChrCount2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad2.vcf", "CHR_COUNTS"),
                0,
                UserException.MalformedFile.class
        );

        executeTest("test bad chr counts #2", spec);
    }

    @Test
    public void testBadID() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "IDS") + " --dbsnp " + b36dbSNP129,
                0,
                UserException.MalformedFile.class
        );

        executeTest("test bad RS ID", spec);
    }

    @Test
    public void testBadAllele() {
        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString("validationExampleBad.vcf", "ALLELES"),
            0,
            UserException.MalformedFile.class
        );

        executeTest("test bad alt allele", spec);
    }

    @Test
    public void testBadAllele2() {
        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString("validationExampleBad3.vcf", "REF"),
            0,
            UserException.MalformedFile.class
        );

        executeTest("test bad ref allele in deletion", spec);
    }

}
