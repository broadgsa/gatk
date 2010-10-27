package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class
        VariantEvalIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf -reportType Grep";

    private static String rootGZ = cmdRoot +
                " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
                " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf.gz" +
                " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf.gz -reportType Grep";

    private static String[] testsEnumerations = {root, rootGZ};

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -o %s",
                    1, Arrays.asList("158ac8e6d32eb2ea1bbeebfa512965de"));
            //executeTestParallel("testSelect1", spec);
            executeTest("testSelect1", spec);
        }
    }

    @Test
    public void testSelect2() {
        String extraArgs = "-L 1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
                1, Arrays.asList("cee96f61ffa1d042fe0c63550c508ec9"));
        //executeTestParallel("testSelect2", spec);
        executeTest("testSelect2", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFiles[] = {"GenotypeConcordanceEval.vcf", "GenotypeConcordanceEval.vcf.gz"};
        for (String vcfFile : vcfFiles) {
            WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -B:eval,VCF " + validationDataLocation + vcfFile + " -B:comp,VCF " + validationDataLocation + "GenotypeConcordanceComp.vcf -noStandard -E GenotypeConcordance -reportType CSV -o %s",
                    1,
                    Arrays.asList("7e9ce1b26cdeaa50705f5de163847638"));
            //executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
            executeTest("testVEGenotypeConcordance" + vcfFile, spec);
        }

    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "ff8c4ba16c7c14b4edbaf440f20641f9");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0 -E MendelianViolationEvaluator", "9b509ce5d31658eb09bb9597799b2908");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();
            for (String tests : testsEnumerations) {
                WalkerTestSpec spec = new WalkerTestSpec( tests + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
                //executeTestParallel("testVESimple", spec);
                executeTest("testVESimple", spec);
            }
        }
    }

    @Test
    public void testVEComplex() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -MVQ 30 -E MendelianViolationEvaluator" +
                " -B:dbsnp_130,dbSNP " + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B:comp_hapmap,VCF " + validationDataLocation + "CEU_hapmap_nogt_23.vcf";


        String matchingMD5 = "1e6d6e152c9a90513dd5b6eaa3729388";
        expectations.put("", matchingMD5);
        expectations.put(" -known comp_hapmap -known dbsnp", matchingMD5);
        expectations.put(" -known comp_hapmap", "d28dd504017f39a91cde8e6f096879d6");
        for (String tests : testsEnumerations) {
            for (Map.Entry<String, String> entry : expectations.entrySet()) {
                String extraArgs2 = entry.getKey();
                String md5 = entry.getValue();

                WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs1 + extraArgs2 + " -o %s",
                        1, // just one output file
                        Arrays.asList(md5));
                //executeTestParallel("testVEComplex", spec);
                executeTest("testVEComplex", spec);
            }
        }
    }

    @Test
    public void testVEGenomicallyAnnotated() {
        String vecmd = "-T VariantEval" +
                       " -R " + b36KGReference +
                       " -L 21" +
                       " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
                       " -E CountFunctionalClasses -noStandard" +
                       " -B:eval,VCF " + validationDataLocation + "test.filtered.maf_annotated.vcf" +
                       " -o %s";
        String md5 = "d41d8cd98f00b204e9800998ecf8427e";

        WalkerTestSpec spec = new WalkerTestSpec(vecmd, 1, Arrays.asList(md5));
        //executeTestParallel("testVEGenomicallyAnnotated", spec);
        executeTest("testVEGenomicallyAnnotated", spec);
    }

    @Test
    public void testVEWriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -NO_HEADER -family NA19238+NA19239=NA19240 -MVQ 30 -E MendelianViolationEvaluator";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs + " -o %s -outputVCF %s -NO_HEADER",
                    2,
                    Arrays.asList("6b97a019402b3984fead9a4e8b7c7c2a", "989bc30dea6c8a4cf771cd1b9fdab488"));
            //executeTestParallel("testVEWriteVCF", spec);
            executeTest("testVEWriteVCF", spec);
        }
    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -E GenotypeConcordance -B:evalYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.ug.very.few.lines.vcf -B:compYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.fake.genotypes.ac.test.vcf -reportType CSV";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("25a681855cb26e7380fbf1a93de0a41f"));
        //executeTestParallel("testACDiscordanceAtAC1EvalAC2Comp",spec);
        executeTest("testACDiscordanceAtAC1EvalAC2Comp",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }
}
