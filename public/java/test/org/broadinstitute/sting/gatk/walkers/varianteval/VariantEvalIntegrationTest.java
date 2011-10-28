package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = validationDataLocation + "VariantEval";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.HG00625.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    @Test
    public void testFunctionClassWithSnpeff() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + validationDataLocation + "snpEff.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-noST",
                                        "-ST FunctionalClass",
                                        "-L " + validationDataLocation + "snpEff.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("d9dcb352c53106f54fcc981f15d35a90")
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
                                Arrays.asList("6a71b17c19f5914c277a99f45f5d9c39")
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
                                Arrays.asList("1fefd6cf9c2554d5f886c3998defd4d0")
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
                Arrays.asList("d470e00a368b5a0468012818994c6a89")
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
                Arrays.asList("12856e52c2682328f91594089328596c")
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
                Arrays.asList("91610b9240f64e0eb03cfd2602cf57af")
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
                Arrays.asList("e40b77e7ed6581328e373a24b93cd170")
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
                Arrays.asList("15beaf3823c131cabc5fb0445239f978")
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
                Arrays.asList("7ddd4ee74938d229ce5cb7b9b9104abe")
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
                Arrays.asList("a90f33906a732ef5eb346e559c96ccc1")
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
                Arrays.asList("2567f90d3d7354850c5a59730ecc6e4f")
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
                Arrays.asList("fa091aa8967893389c51102fd9f0bebb")
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
                1, Arrays.asList("f70997b6a3e7fdc89d11e1d61a2463d4"));
        executeTestParallel("testSelect1", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFile = "GenotypeConcordanceEval.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -ST CpG --eval:VCF3 " + validationDataLocation + vcfFile + " --comp:VCF3 " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                1,
                Arrays.asList("96f27163f16bb945f19c6623cd6db34e"));
        executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -ST CpG -EV GenotypeConcordance --eval:evalYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.very.few.lines.vcf --comp:compYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("407682de41dcf139ea635e9cda21b912"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("009ecc8376a20dce81ff5299ef6bfecb"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("424c9d438b1faa59b2c29413ba32f37b"));
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
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("18fa0b89ebfff51141975d7e4ce7a159"));
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

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("0b81d97f843ec4a1a4222d1f9949bfca"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults() {
        String md5 = "b0565ac61b2860248e4abd478a177b5e";

        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-sn HG00625",
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
                                Arrays.asList("da65fc8f0d0eeaf0a0b06a07f444bb8e")
                              );
        executeTest("testAlleleCountStrat", spec);
    }
}
