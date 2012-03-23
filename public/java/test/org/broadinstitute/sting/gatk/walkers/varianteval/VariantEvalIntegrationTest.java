/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = validationDataLocation + "VariantEval";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.NA12045.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    @Test
    public void testFunctionClassWithSnpeff() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + validationDataLocation + "snpEff2.0.5.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-noST",
                                        "-ST FunctionalClass",
                                        "-L " + validationDataLocation + "snpEff2.0.5.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("add8b2213c091a41f5d7a2c8dd68c03a")
                              );
        executeTest("testFunctionClassWithSnpeff", spec);
    }

    @Test
    public void testStratifySamplesAndExcludeMonomorphicSites() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + variantEvalTestDataRoot + "/CEU.trio.callsForVE.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-ST Sample",
                                        "-L " + variantEvalTestDataRoot + "/CEU.trio.callsForVE.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("621a712deb01e7fc7e5a13d3627b11ba")
                              );
        executeTest("testStratifySamplesAndExcludeMonomorphicSites", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + fundamentalTestVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-L " + fundamentalTestVCF,
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("7a726ecbedd722fa7cd4de3e023b7a82")
                              );
        executeTest("testFundamentalsCountVariantsSNPsandIndels", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("95bb4a4267a8f29dd7a8169561499f20")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNovelty", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-ST Filter",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("9b51029083495935823fb0447a2857b9")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNoveltyAndFilter", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST CpG",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("318b5fbbc61e2fc11d49369359812edd")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithCpG", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST FunctionalClass",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("74c02df2ef69dda231a2aec2a948747b")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithFunctionalClass", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Degeneracy",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("2d97b1fe15e532e89803ba7ba347ff20")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithDegeneracy", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Sample",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("474cbc231ddbc4ba299ffe61a17405b6")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithSample", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthSelect",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("2cc9bc4bbe8b4edb6dc27642ec41f66e")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithJexlExpression", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthLt20",
                        "-select 'DP > 20'",
                        "-selectName DepthGt20",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("00c94cf3e14bc2855d39bbefa27f9bb2")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithMultipleJexlExpressions", spec);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("a0c0d4805db1245aa30a306aa506096f")
        );
        executeTest("testFundamentalsCountVariantsNoCompRod", spec);
    }

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        String tests = cmdRoot +
                " --dbsnp " + b36dbSNP129 +
                " --eval " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
                " --comp:comp_genotypes,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.head.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -ST CpG -o %s",
                1, Arrays.asList("2192418a70a8e018a1675d4f425155f3"));
        executeTestParallel("testSelect1", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFile = "GenotypeConcordanceEval.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -ST CpG --eval:VCF3 " + validationDataLocation + vcfFile + " --comp:VCF3 " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                1,
                Arrays.asList("70da6a0f91a9f1052d68fc360cc99aed"));
        executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
    }

    @Test
    public void testVEMendelianViolationEvaluator() {
        String vcfFile = "/MendelianViolationEval.vcf";
        String pedFile = "/MendelianViolationEval.ped";

        WalkerTestSpec spec = new WalkerTestSpec("-T VariantEval -R "+b37KGReference+" --eval " + variantEvalTestDataRoot + vcfFile + " -ped "+ variantEvalTestDataRoot + pedFile +" -noEV -EV MendelianViolationEvaluator -L 1:10109-10315 -o %s -mvq 0 -noST",
                1,
                Arrays.asList("03581adcb4f2f7960662fc7ffd910f43"));
        executeTestParallel("testVEMendelianViolationEvaluator" + vcfFile, spec);
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -ST CpG -EV GenotypeConcordance --eval:evalYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.very.few.lines.vcf --comp:compYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("2282523336c24d434d1cc0eb1697b4f9"));
        executeTestParallel("testCompVsEvalAC",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test(enabled = false) // no longer supported in the GATK
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" --eval " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -ST CpG -tf " + testDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("6af2b9959aa1778a5b712536de453952"));
        executeTestParallel("testTranches",spec);
    }

    @Test
    public void testCompOverlap() {
        String extraArgs = "-T VariantEval -R " + b37KGReference + " -L " + validationDataLocation + "VariantEval/pacbio.hg19.intervals --comp:comphapmap " + comparisonDataLocation + "Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf --eval " + validationDataLocation + "VariantEval/pacbio.ts.recalibrated.vcf -noEV -EV CompOverlap -sn NA12878 -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("00241ce70476187a2f910606b9242697"));
        executeTestParallel("testCompOverlap",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " +
                           b37KGReference +
                           " -L 20" +
                           " --dbsnp " + b37dbSNP132 +
                           " --eval:evalBI " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                           " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("ec321fcc424fbad74a4a74e739173d03"));
        executeTestParallel("testEvalTrackWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " + b37KGReference +
                " -L 20" +
                " --dbsnp " + b37dbSNP132 +
                " --eval:evalBI " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                " --eval:evalBC " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bc.sites.vcf" +
                " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("ccaea6245086552cd63f828eabddfaf3"));
        executeTestParallel("testMultipleEvalTracksWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleCompTracks() {
        String dbsnp = GATKDataLocation + "dbsnp_132_b37.vcf";

        String extraArgs =  "-T VariantEval" +
                           " -R " + b37KGReference +
                           " --comp " + validationDataLocation + "/VariantEval/ALL.phase1.chr20.broad.snps.genotypes.subset.vcf" +
                           " --eval " + validationDataLocation + "/VariantEval/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf" +
                           " --dbsnp " + dbsnp +
                           " -L 20:10000000-10100000" +
                           " -noST -noEV -ST Novelty -EV CompOverlap" +
                           " -o %s";

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("15f6a6ba4f7fed49c617589ce9fdcbc5"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults1() {
        String md5 = "bcf55537db0762b8fd68f7f02439c475";

        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-sn NA12045",
                        "-noST",
                        "-L " + fundamentalTestSNPsVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList(md5)
        );
        executeTestParallel("testPerSampleAndSubsettedSampleHaveSameResults-subset", spec);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsOneSampleVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-noST",
                        "-L " + fundamentalTestSNPsOneSampleVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList(md5)
        );
        executeTestParallel("testPerSampleAndSubsettedSampleHaveSameResults-onesample", spec2);
    }


    @Test
    public void testAlleleCountStrat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + fundamentalTestSNPsVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-ST AlleleCount",
                                        "-L " + fundamentalTestSNPsVCF,
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("9954c769ef37c47d3b61481ab0807be0")
                              );
        executeTest("testAlleleCountStrat", spec);
    }

    @Test
    public void testIntervalStrat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-eval " + testDir + "/withSymbolic.b37.vcf",
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-stratIntervals " + testDir + "/overlapTest.bed",
                                        "-ST IntervalStratification",
                                        "-L 20",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("c0d69ce7647a575d166d8bab5aa16299")
                              );
        executeTest("testIntervalStrat", spec);
    }

    @Test
    public void testModernVCFWithLargeIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-eval " + validationDataLocation + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.vcf",
                                        "-L 20",
                                        "-D " + b37dbSNP132,
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("9a8ffb506118c1bde6f7bfadc4fb6f10")
                              );
        executeTest("testModernVCFWithLargeIndels", spec);
    }
}
