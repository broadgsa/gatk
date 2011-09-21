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
import org.testng.annotations.Test;

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
                 baseTestString(" -priority v1 -V:v1 " + validationDataLocation + file + args),
                 1,
                 Arrays.asList(md5));
         executeTest("testInOut1--" + file, spec);
    }

    public void combine2(String file1, String file2, String args, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1,v2 -V:v1 " + validationDataLocation + file1 + " -V:v2 "+ validationDataLocation + file2 + args),
                 1,
                 Arrays.asList(md5));
         executeTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }

    public void combineSites(String args, String md5) {
        String file1 = "1000G_omni2.5.b37.sites.vcf";
        String file2 = "hapmap_3.3.b37.sites.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants -NO_HEADER -o %s -R " + b37KGReference
                        + " -L 1:1-10,000,000 -V:omni " + validationDataLocation + file1
                        + " -V:hm3 " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    public void combinePLs(String file1, String file2, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants -NO_HEADER -o %s -R " + b36KGReference + " -priority v1,v2 -V:v1 " + validationDataLocation + file1 + " -V:v2 " + validationDataLocation + file2,
                 1,
                 Arrays.asList(md5));
         executeTest("combine PLs 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }

    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "c608b9fc1e36dba6cebb4f259883f9f0"); }
    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "20caad94411d6ab48153b214de916df8", " -setKey foo"); }
    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "004f3065cb1bc2ce2f9afd695caf0b48", " -setKey null"); }
    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "c9c901ff9ef2a982624b203a8086dff0"); } // official project VCF files in tabix format

    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "7593be578d4274d672fc22fced38012b"); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "1cd467863c4e948fadd970681552d57e"); }

    @Test public void combineWithPLs() { combinePLs("combine.3.vcf", "combine.4.vcf", "0f873fed02aa99db5b140bcd6282c10a"); }

    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "1d5a021387a8a86554db45a29f66140f"); } // official project VCF files in tabix format
    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "20163d60f18a46496f6da744ab5cc0f9"); } // official project VCF files in tabix format
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "d76cd5b3ced7745d42fe0af39ce0b32e"); }

    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "e144b6283765494bfe8189ac59965083"); }

    @Test public void uniqueSNPs() { combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "", "89f55abea8f59e39d1effb908440548c"); }

    @Test public void omniHM3Union() { combineSites(" -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED", "c6adeda751cb2a08690dd9202356629f"); }
    @Test public void omniHM3Intersect() { combineSites(" -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED", "3a08fd5ee18993dfc8882156ccf5d2e9"); }

    @Test public void threeWayWithRefs() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:NA19240_BGI "+validationDataLocation+"NA19240.BGI.RG.vcf" +
                        " -V:NA19240_ILLUMINA "+validationDataLocation+"NA19240.ILLUMINA.RG.vcf" +
                        " -V:NA19240_WUGSC "+validationDataLocation+"NA19240.WUGSC.RG.vcf" +
                        " -V:denovoInfo "+validationDataLocation+"yri_merged_validation_data_240610.annotated.b36.vcf" +
                        " -setKey centerSet" +
                        " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
                        " -priority NA19240_BGI,NA19240_ILLUMINA,NA19240_WUGSC,denovoInfo" +
                        " -genotypeMergeOptions UNIQUIFY -L 1"),
                1,
                Arrays.asList("1de95f91ca15d2a8856de35dee0ce33e"));
        executeTest("threeWayWithRefs", spec);
    }

    // complex examples with filtering, indels, and multiple alleles
    public void combineComplexSites(String args, String md5) {
        String file1 = "combine.1.vcf";
        String file2 = "combine.2.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants -NO_HEADER -o %s -R " + b37KGReference
                        + " -V:one " + validationDataLocation + file1
                        + " -V:two " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineComplexSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    @Test public void complexTestFull() { combineComplexSites("", "b5a53ee92bdaacd2bb3327e9004ae058"); }
    @Test public void complexTestMinimal() { combineComplexSites(" -minimalVCF", "df96cb3beb2dbb5e02f80abec7d3571e"); }
    @Test public void complexTestSitesOnly() { combineComplexSites(" -sites_only", "f704caeaaaed6711943014b847fe381a"); }
    @Test public void complexTestSitesOnlyMinimal() { combineComplexSites(" -sites_only -minimalVCF", "f704caeaaaed6711943014b847fe381a"); }

    @Test
    public void combineDBSNPDuplicateSites() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants -NO_HEADER -L 1:902000-903000 -o %s -R " + b37KGReference + " -V:v1 " + b37dbSNP132,
                 1,
                 Arrays.asList(""));
         executeTest("combineDBSNPDuplicateSites:", spec);
    }
}