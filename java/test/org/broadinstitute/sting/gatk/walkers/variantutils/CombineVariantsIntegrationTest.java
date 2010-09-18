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
        return "-T CombineVariants -NO_HEADER -L 1:1-50,000,000 -o %s -R " + b36KGReference + args;
    }

    public void test1InOut(String file, String md5) {
        test1InOut(file, md5, "");
    }

    public void test1InOut(String file, String md5, String args) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1 -B:v1,VCF " + validationDataLocation + file + args),
                 1,
                 Arrays.asList(md5));
         executeTest("testInOut1--" + file, spec);
    }

    public void combine2(String file1, String file2, String args, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1,v2 -B:v1,VCF " + validationDataLocation + file1 + " -B:v2,VCF " + validationDataLocation + file2 + args),
                 1,
                 Arrays.asList(md5));
         executeTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }


    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "def6365e318cb506aa69ee8ec6899c4e"); }
    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "6d1d9fffcc10d4f9ab086760603b577e", " -setKey foo"); }
    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "e1d43e0c4df65df5e9761916b12f0b5b", " -setKey null"); }
    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "79e7e474b0a2bb82930201f143328d5d"); } // official project VCF files in tabix format

    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "d93ee8b9f36eedc80a3afa548bffc888"); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "4d16ffe9cf3beff94ac05ea49f85049a"); }

    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "78f2a2c09df8fa51ccfe5b3706a5e313"); } // official project VCF files in tabix format
    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "460bdbeaf7e9641395ac2ce6e1afc106"); } // official project VCF files in tabix format
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "d445583e7f3a8bc98004e8e44878fce6"); }

    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "e64284f598fc4f473ee5554282d84b0c"); }

    @Test public void uniqueSNPs() { combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "", "1f3d118dcee9b711714856639930694d"); }

    @Test public void threeWayWithRefs() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -B:NA19240_BGI,VCF "+validationDataLocation+"NA19240.BGI.RG.vcf" +
                        " -B:NA19240_ILLUMINA,VCF "+validationDataLocation+"NA19240.ILLUMINA.RG.vcf" +
                        " -B:NA19240_WUGSC,VCF "+validationDataLocation+"NA19240.WUGSC.RG.vcf" +
                        " -B:denovoInfo,VCF "+validationDataLocation+"yri_merged_validation_data_240610.annotated.b36.vcf" +
                        " -setKey centerSet" +
                        " -variantMergeOptions UNION" +
                        " -priority NA19240_BGI,NA19240_ILLUMINA,NA19240_WUGSC,denovoInfo" +
                        " -genotypeMergeOptions UNIQUIFY -L 1"),
                1,
                Arrays.asList("81f53a9f2a6d410c1236266e2516bf6f"));
        executeTest("threeWayWithRefs", spec);
        
    }
}
