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
                                        "--eval " + validationDataLocation + "snpEff2.0.4.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-noST",
                                        "-ST FunctionalClass",
                                        "-L " + validationDataLocation + "snpEff2.0.4.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("f909fd8374f663e983b9b3fda4cf5cf1")
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
                                Arrays.asList("081fcaa532c7ba8f23da739389e6f7c3")
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
                                Arrays.asList("b3852f84d07c270b8a12874083c3e31b")
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
                Arrays.asList("cf70468b5ebaec408419da69b0a7fcb9")
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
                Arrays.asList("5e3b8b85acfc41365c8208c23abf746b")
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
                Arrays.asList("ccdbc50d30ece6d0d3b199c397f03ed3")
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
                Arrays.asList("95c690d5af8ed51573eb2f0503dcd9c2")
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
                Arrays.asList("8e8547eb38b34bec0095b0500fd9641d")
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
                Arrays.asList("158a4651a656aea7f84c79548f6fe519")
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
                Arrays.asList("76c8a0b28d2993644120f7afa5833ab2")
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
                Arrays.asList("34682193f458b93b39efac00b4fc6723")
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
                Arrays.asList("52f6655f1532bcea24b402010d93ce73")
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
                1, Arrays.asList("4f60acc8a4b21c4b4ccb51ad9071449c"));
        executeTestParallel("testSelect1", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFile = "GenotypeConcordanceEval.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -ST CpG --eval:VCF3 " + validationDataLocation + vcfFile + " --comp:VCF3 " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                1,
                Arrays.asList("9a56c20a7b9a554a7b530f2cb1dd776d"));
        executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -ST CpG -EV GenotypeConcordance --eval:evalYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.very.few.lines.vcf --comp:compYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("fa13eb59892892c07711c6ffe31bf870"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("9002023b8aa8fc2c9aac58b8a79bca1e"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("190e1a171132832bf92fbca56a9c40bb"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("08586d443fdcf3b7f63b8f9e3a943c62"));
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

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("61052c19211e7eb61fbbb62db5e40b56"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults1() {
        String md5 = "0edded1cd578db62fa296c99c34a909d";

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
                                Arrays.asList("ee22604616b3e9fc48a6dcbbf73a056d")
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
                                Arrays.asList("240369cd651c77e05e8a6659f4a6237e")
                              );
        executeTest("testIntervalStrat", spec);
    }
}
