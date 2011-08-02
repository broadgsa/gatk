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

    public void test1InOut(String file, String md5, boolean vcf3) {
        test1InOut(file, md5, "", vcf3);
    }

    public void test1InOut(String file, String md5, String args, boolean vcf3) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1 -B:v1,VCF" + (vcf3 ? "3 " : " ") + validationDataLocation + file + args),
                 1,
                 Arrays.asList(md5));
         executeTest("testInOut1--" + file, spec);
    }

    public void combine2(String file1, String file2, String args, String md5, boolean vcf3) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1,v2 -B:v1,VCF" + (vcf3 ? "3 " : " ") + validationDataLocation + file1 + " -B:v2,VCF" + (vcf3 ? "3 " : " ") + validationDataLocation + file2 + args),
                 1,
                 Arrays.asList(md5));
         executeTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec);
    }

    public void combineSites(String args, String md5) {
        String file1 = "1000G_omni2.5.b37.sites.vcf";
        String file2 = "hapmap_3.3.b37.sites.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants -NO_HEADER -o %s -R " + b37KGReference
                        + " -L 1:1-10,000,000 -B:omni,VCF " + validationDataLocation + file1
                        + " -B:hm3,VCF " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "c608b9fc1e36dba6cebb4f259883f9f0", true); }
    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "20caad94411d6ab48153b214de916df8", " -setKey foo", true); }
    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "004f3065cb1bc2ce2f9afd695caf0b48", " -setKey null", true); }
    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "c9c901ff9ef2a982624b203a8086dff0", false); } // official project VCF files in tabix format

    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "7593be578d4274d672fc22fced38012b", false); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "1cd467863c4e948fadd970681552d57e", false); }

    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "1d5a021387a8a86554db45a29f66140f", false); } // official project VCF files in tabix format
    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "20163d60f18a46496f6da744ab5cc0f9", false); } // official project VCF files in tabix format
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "f1cf095c2fe9641b7ca1f8ee2c46fd4a", false); }

    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "e144b6283765494bfe8189ac59965083", false); }

    @Test public void uniqueSNPs() { combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "", "89f55abea8f59e39d1effb908440548c", true); }

    @Test public void omniHM3Union() { combineSites(" -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED", "4836086891f6cbdd40eebef3076d215a"); }
    @Test public void omniHM3Intersect() { combineSites(" -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED", "6a34b5d743efda8b2f3b639f3a2f5de8"); }

    @Test public void threeWayWithRefs() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -B:NA19240_BGI,VCF "+validationDataLocation+"NA19240.BGI.RG.vcf" +
                        " -B:NA19240_ILLUMINA,VCF "+validationDataLocation+"NA19240.ILLUMINA.RG.vcf" +
                        " -B:NA19240_WUGSC,VCF "+validationDataLocation+"NA19240.WUGSC.RG.vcf" +
                        " -B:denovoInfo,VCF "+validationDataLocation+"yri_merged_validation_data_240610.annotated.b36.vcf" +
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
                        + " -B:one,VCF " + validationDataLocation + file1
                        + " -B:two,VCF " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        executeTest("combineComplexSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec);
    }

    @Test public void complexTestFull() { combineComplexSites("", "b5a53ee92bdaacd2bb3327e9004ae058"); }
    @Test public void complexTestMinimal() { combineComplexSites(" -minimalVCF", "df96cb3beb2dbb5e02f80abec7d3571e"); }
    @Test public void complexTestSitesOnly() { combineComplexSites(" -sites_only", "f72a178137e25dbe0b931934cdc0079d"); }
    @Test public void complexTestSitesOnlyMinimal() { combineComplexSites(" -sites_only -minimalVCF", "f704caeaaaed6711943014b847fe381a"); }
}