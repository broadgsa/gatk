package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = validationDataLocation + "/VariantEval";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.HG00625.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-B:dbsnp " + b37dbSNP132,
                                        "-B:eval " + fundamentalTestVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-BTI eval",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("bced1842c78fbabb089dd12b7087050d")
                              );
        executeTest("testFundamentalsCountVariantsSNPsandIndels", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("06510bd37ffaa39e817ca0dcaf8f8ac2")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNovelty", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-ST Filter",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("19c5b1b6396921c5b1059a2849ae4fcc")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNoveltyAndFilter", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST CpG",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("a71f8d81cf166cd97ac628092650964a")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithCpG", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST FunctionalClass",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("4dabe0658232f6174188515db6dfe112")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithFunctionalClass", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Degeneracy",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("3340587f10ceff83e5567ddfd1a9a60e")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithDegeneracy", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Sample",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("c730c7ee31c8138cef6efd8dd04fbbfc")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithSample", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthSelect",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("2559ca8f454b03e81561f6947f79df18")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithJexlExpression", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthLt20",
                        "-select 'DP > 20'",
                        "-selectName DepthGt20",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("23aa5f97641d2fd033095f21c51d2f37")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithMultipleJexlExpressions", spec);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("a69dd3f06903b3f374c6d6f010c653e0")
        );
        executeTest("testFundamentalsCountVariantsNoCompRod", spec);
    }

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        String tests = cmdRoot +
                " -B:dbsnp " + b36dbSNP129 +
                " -B:eval,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
                " -B:comp_genotypes,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.head.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -ST CpG -o %s",
                1, Arrays.asList("14054badcd89b24c2375e1d09918f681"));
        executeTestParallel("testSelect1", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFile = "GenotypeConcordanceEval.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -ST CpG -B:eval,VCF3 " + validationDataLocation + vcfFile + " -B:comp,VCF3 " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                1,
                Arrays.asList("96f27163f16bb945f19c6623cd6db34e"));
        executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -ST CpG -EV GenotypeConcordance -B:evalYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.very.few.lines.vcf -B:compYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("d1932be3748fcf6da77dc51aec323710"));
        executeTestParallel("testCompVsEvalAC",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" -B:eval,vcf " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -ST CpG -tf " + testDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("984df6e94a546294fc7e0846cbac2dfe"));
        executeTestParallel("testTranches",spec);
    }

    @Test
    public void testCompOverlap() {
        String extraArgs = "-T VariantEval -R " + b37KGReference + " -L " + validationDataLocation + "VariantEval/pacbio.hg19.intervals -B:comphapmap,vcf " + comparisonDataLocation + "Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf -B:eval,vcf " + validationDataLocation + "VariantEval/pacbio.ts.recalibrated.vcf -noEV -EV CompOverlap -sn NA12878 -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("462d4784dd55294ef9d5118217b157a5"));
        executeTestParallel("testCompOverlap",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " +
                           b37KGReference +
                           " -L 20" +
                           " -B:dbsnp " + b37dbSNP132 +
                           " -B:evalBI " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                           " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("0897dfba2f4a245faddce38000555cce"));
        executeTestParallel("testEvalTrackWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " + b37KGReference +
                " -L 20" +
                " -B:dbsnp " + b37dbSNP132 +
                " -B:evalBI " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                " -B:evalBC " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bc.sites.vcf" +
                " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("ead3602e14ec2944b5d9e4dacc08c819"));
        executeTestParallel("testMultipleEvalTracksWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleCompTracks() {
        String dbsnp = GATKDataLocation + "dbsnp_132_b37.vcf";

        String extraArgs =  "-T VariantEval" +
                           " -R " + b37KGReference +
                           " -B:comp " + validationDataLocation + "/VariantEval/ALL.phase1.chr20.broad.snps.genotypes.subset.vcf" +
                           " -B:eval " + validationDataLocation + "/VariantEval/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf" +
                           " -B:dbsnp " + dbsnp +
                           " -L 20:10000000-10100000" +
                           " -noST -noEV -ST Novelty -EV CompOverlap" +
                           " -o %s";

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("167a347ce0729d1bc3d4fd5069ebd674"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults() {
        String md5 = "40471a84b501eb440ee2d42e3081f228";

        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestSNPsVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-sn HG00625",
                        "-noST",
                        "-BTI eval",
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
                        "-B:dbsnp " + b37dbSNP132,
                        "-B:eval " + fundamentalTestSNPsOneSampleVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-noST",
                        "-BTI eval",
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
                                        "-B:dbsnp " + b37dbSNP132,
                                        "-B:eval " + fundamentalTestSNPsVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-ST AlleleCount",
                                        "-BTI eval",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("44464fe7c89a56cf128a932ef640f7da")
                              );
        executeTest("testAlleleCountStrat", spec);
    }
}
