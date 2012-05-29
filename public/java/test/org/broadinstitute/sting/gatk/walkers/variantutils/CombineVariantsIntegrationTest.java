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
        return "-T CombineVariants --no_cmdline_in_header -L 1:1-50,000,000 -o %s -R " + b36KGReference + args;
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
                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
                        + " -L 1:1-10,000,000 -V:omni " + validationDataLocation + file1
                        + " -V:hm3 " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    public void combinePLs(String file1, String file2, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants --no_cmdline_in_header -o %s -R " + b36KGReference + " -priority v1,v2 -V:v1 " + validationDataLocation + file1 + " -V:v2 " + validationDataLocation + file2,
                 1,
                 Arrays.asList(md5));
         executeTest("combine PLs 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }

    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "c1e82f0842ca721d10f21604f26a5248"); }
    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "b2fcf3983cc9e667b9bbed8372080776", " -setKey foo"); }
    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "98c0cbb94e5debf7545a656665a1b659", " -setKey null"); }
    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "10170f9e72cc831a5820bd03e70fe46a"); } // official project VCF files in tabix format

    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "074e909f80ffcc9fddc3fac89ea36bef"); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "f26980af214011c0452b8ce843f3063b"); }

    @Test public void combineWithPLs() { combinePLs("combine.3.vcf", "combine.4.vcf", "5bc1de1197506aced0f9e7a08b572c44"); }

    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "06e86711bcf0efc0f0c4a378f6147cf6"); } // official project VCF files in tabix format
    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "03103f6b39e9fb7a396df0013f01fae6"); } // official project VCF files in tabix format
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "12fc1b8145f7884762f0c2cbbd319ae1"); }

    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "7e2dba80ba38b2a86713f635d630eb59"); }

    @Test public void uniqueSNPs() { combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "", "3950392e1b8b53ae363e705185ad1da9"); }

    @Test public void omniHM3Union() { combineSites(" -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED", "5c60eb8d5d4b957a0cf52ca008f021ba"); }
    @Test public void omniHM3Intersect() { combineSites(" -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED", "774b43e69cc7ec93090b4f6e9f4a1079"); }

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
                Arrays.asList("948291bbf47d1cec692d0fe4358ff92c"));
        executeTest("threeWayWithRefs", spec);
    }

    // complex examples with filtering, indels, and multiple alleles
    public void combineComplexSites(String args, String md5) {
        String file1 = "combine.1.vcf";
        String file2 = "combine.2.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
                        + " -V:one " + validationDataLocation + file1
                        + " -V:two " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineComplexSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    @Test public void complexTestFull() { combineComplexSites("", "dd805f6edfc3cf724512dfbbe8df5183"); }
    @Test public void complexTestMinimal() { combineComplexSites(" -minimalVCF", "14a205edb022f79abf1863588cfee56b"); }
    @Test public void complexTestSitesOnly() { combineComplexSites(" -sites_only", "e118d04d1d47c02ad38c046561a9f616"); }
    @Test public void complexTestSitesOnlyMinimal() { combineComplexSites(" -sites_only -minimalVCF", "e118d04d1d47c02ad38c046561a9f616"); }

    @Test
    public void combineDBSNPDuplicateSites() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants --no_cmdline_in_header -L 1:902000-903000 -o %s -R " + b37KGReference + " -V:v1 " + b37dbSNP132,
                 1,
                 Arrays.asList("a838dc241cf357466cd4331fd298c73a"));
         executeTest("combineDBSNPDuplicateSites:", spec);
    }
}