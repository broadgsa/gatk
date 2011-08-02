package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = validationDataLocation + "/VariantEval";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = variantEvalTestDataRoot + "/" + "FundamentalsTest.annotated.db.subset.final.HG00625.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -B:comp_genotypes,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.head.vcf";

    private static String rootGZ = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf.gz" +
            " -B:comp_genotypes,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.head.vcf.gz";

    // TODO -- I can't seem to reindex this VCF using Tabix without it causing failures.  Looking into it.  [EB]
    // private static String[] testsEnumerations = {root, rootGZ};
    private static String[] testsEnumerations = {root};

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                                        "-B:eval,VCF " + fundamentalTestVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-BTI eval",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("48b8417c1f8bd74ff7b9808580abd2a2")
                              );
        executeTest("testFundamentalsCountVariantsSNPsandIndels", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("86d45ecefdf5849c55b3ca8f82a3d525")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNovelty", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-ST Filter",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("3d18901ec1766aa2e748eac913f5ddcd")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNoveltyAndFilter", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST CpG",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("677fe398643e62a10d6739d36a720a12")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithCpG", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST FunctionalClass",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("5fb44fd7cb00941c986a9941e43e44cd")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithFunctionalClass", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Degeneracy",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("daaca7ef3b7313e5af217cbc6f37c9e2")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithDegeneracy", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Sample",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("97c466f8ffd0fcf2c30ef08669d213d9")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithSample", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
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
                Arrays.asList("df8cdfcf3d0c2fc795812c6eae6a76f8")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithJexlExpression", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestVCF,
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
                Arrays.asList("c7aed12265e2b2311d17a0cc8a29f6aa")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithMultipleJexlExpressions", spec);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:eval,VCF " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-BTI eval",
                        "-o %s"
                ),
                1,
                Arrays.asList("d44c8f44384189a09eea85a8e89d7299")
        );
        executeTest("testFundamentalsCountVariantsNoCompRod", spec);
    }

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -ST CpG -o %s",
                    1, Arrays.asList("96860dedea0fa6b46c07f46b847fea42"));
            executeTestParallel("testSelect1", spec);
        }
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFile = "GenotypeConcordanceEval.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -ST CpG -B:eval,VCF3 " + validationDataLocation + vcfFile + " -B:comp,VCF3 " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                1,
                Arrays.asList("e4c981f7f5d78680c71310fc9be9a1c1"));
        executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -ST CpG -EV GenotypeConcordance -B:evalYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.ug.very.few.lines.vcf -B:compYRI,VCF3 " + validationDataLocation + "yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("162daa5039e1965eb2423a8589339a69"));
        executeTestParallel("testCompVsEvalAC",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" -B:eval,vcf " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -ST CpG -tf " + testDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("90cd98044e754b80034a9f4e6d2c55b9"));
        executeTestParallel("testTranches",spec);
    }

    @Test
    public void testCompOverlap() {
        String extraArgs = "-T VariantEval -R " + b37KGReference + " -L " + validationDataLocation + "VariantEval/pacbio.hg19.intervals -B:comphapmap,vcf " + comparisonDataLocation + "Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf -B:eval,vcf " + validationDataLocation + "VariantEval/pacbio.ts.recalibrated.vcf -noEV -EV CompOverlap -sn NA12878 -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("70aa420929de7f888a6f48c2d01bbcda"));
        executeTestParallel("testCompOverlap",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() {
        String dbsnp = GATKDataLocation + "dbsnp_129_b37.rod";

        String extraArgs = "-T VariantEval -R " +
                           b37KGReference +
                           " -L 20" +
                           " -D " + dbsnp +
                           " -B:evalBI,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                           " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("58fdc6c42fade3007537bb99fb3ce738"));
        executeTestParallel("testEvalTrackWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() {
        String dbsnp = GATKDataLocation + "dbsnp_129_b37.rod";

        String extraArgs = "-T VariantEval -R " + b37KGReference +
                " -L 20" +
                " -D " + dbsnp +
                " -B:evalBI,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bi.sites.vcf" +
                " -B:evalBC,VCF " + validationDataLocation + "VariantEval/ALL.20100201.chr20.bc.sites.vcf" +
                " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("34df2815d27e5e62f1694731a7e7953c"));
        executeTestParallel("testMultipleEvalTracksWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleCompTracks() {
        String dbsnp = GATKDataLocation + "dbsnp_132_b37.vcf";

        String extraArgs =  "-T VariantEval" +
                           " -R " + b37KGReference +
                           " -B:comp,VCF " + validationDataLocation + "/VariantEval/ALL.phase1.chr20.broad.snps.genotypes.subset.vcf" +
                           " -B:eval,VCF " + validationDataLocation + "/VariantEval/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf" +
                           " -B:dbsnp,VCF " + dbsnp +
                           " -L 20:10000000-10100000" +
                           " -noST -noEV -ST Novelty -EV CompOverlap" +
                           " -o %s";

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("20332902ae36a84b2fd80405410815f1"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults() {
        String md5 = "9d61f6e2c8592dcf616712a2c587b2af";

        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestSNPsVCF,
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
                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                        "-B:eval,VCF " + fundamentalTestSNPsOneSampleVCF,
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
                                        "-B:dbsnp,VCF " + GATKDataLocation + "dbsnp_132_b37.vcf",
                                        "-B:eval,VCF " + fundamentalTestSNPsVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-ST AlleleCount",
                                        "-BTI eval",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("bf324e4c87fe0d21170fcd2a67a20371")
                              );
        executeTest("testAlleleCountStrat", spec);
    }
}
