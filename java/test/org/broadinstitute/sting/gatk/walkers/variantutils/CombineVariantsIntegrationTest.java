/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Tests CombineVariants
 */
public class CombineVariantsIntegrationTest extends WalkerTest {

    public static String baseTestString(String args) {
        return "-T CombineVariants -L 1:1-50,000,000 -o %s -R " + oneKGLocation + "reference/human_b36_both.fasta" + args;
    }

    // todo -- merge with one file should be identical to input



    public void test1InOut(String file, String type, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1 -B v1," + type + "," + validationDataLocation + file),
                 1,
                 Arrays.asList(md5));
         executeTest("testInOut1--" + file, spec);
    }

    public void combine2(String file1, String type1, String file2, String type2, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1,v2 -B v1," + type1 + "," + validationDataLocation + file1 + " -B v2," + type2 + "," + validationDataLocation + file2),
                 1,
                 Arrays.asList(md5));
         executeTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }


    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "VCF", ""); }
    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "VCF4", ""); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "VCF4", ""); }

    @Test public void combineSNPsAndIndels() { combine2("pilot2.snps.vcf4.genotypes.vcf", "VCF", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "VCF4", ""); }
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "VCF4", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "VCF4", ""); }

    @Test public void uniqueSNPs() { combine2("pilot2.snps.vcf4.genotypes.vcf", "VCF", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "VCF", ""); }
}
