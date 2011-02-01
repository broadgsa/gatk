package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf";

    private static String rootGZ = cmdRoot +
                " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
                " -B:eval,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf.gz" +
                " -B:comp_genotypes,VCF " + validationDataLocation + "yri.trio.gatk.ug.head.vcf.gz";

    private static String[] testsEnumerations = {root, rootGZ};

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -o %s",
                    1, Arrays.asList("8ddc1c4a86cb3f4c22346497785b23e3"));
            //executeTestParallel("testSelect1", spec);
            executeTest("testSelect1", spec);
        }
    }

//    @Test
//    public void testSelect2() {
//        String extraArgs = "-L 1:1-10,000,000";
//        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
//                1, Arrays.asList(""));
//        //executeTestParallel("testSelect2", spec);
//        executeTest("testSelect2", spec);
//    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFiles[] = {"GenotypeConcordanceEval.vcf", "GenotypeConcordanceEval.vcf.gz"};
        for (String vcfFile : vcfFiles) {
            WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -B:eval,VCF " + validationDataLocation + vcfFile + " -B:comp,VCF " + validationDataLocation + "GenotypeConcordanceComp.vcf -noEV -EV GenotypeConcordance -o %s",
                    1,
                    Arrays.asList("7a6754176b573d14b6be7c808a04929d"));
            //executeTestParallel("testVEGenotypeConcordance" + vcfFile, spec);
            executeTest("testVEGenotypeConcordance" + vcfFile, spec);
        }

    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "ffd1abed44faf1590d9026e478b2f8ee");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -mvq 0 -EV MendelianViolationEvaluator", "32c2411fbf58bae5750c8229d15b98eb");

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
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -mvq 30 -EV MendelianViolationEvaluator" +
                " -B:dbsnp_130,dbSNP " + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B:comp_hapmap,VCF " + validationDataLocation + "CEU_hapmap_nogt_23.vcf";


        String matchingMD5 = "6c2fa6573cc57ef8795e9cce2b654d0b";
        expectations.put("", matchingMD5);
        expectations.put(" -knownName comp_hapmap -knownName dbsnp", matchingMD5);
        expectations.put(" -knownName comp_hapmap", "6c2fa6573cc57ef8795e9cce2b654d0b");
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

//    @Test
//    public void testVEGenomicallyAnnotated() {
//        String vecmd = "-T VariantEval" +
//                       " -R " + b36KGReference +
//                       " -L 21" +
//                       " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
//                       " -EV CountFunctionalClasses -noEV" +
//                       " -B:eval,VCF " + validationDataLocation + "test.filtered.maf_annotated.vcf" +
//                       " -o %s";
//        String md5 = "";
//
//        WalkerTestSpec spec = new WalkerTestSpec(vecmd, 1, Arrays.asList(md5));
//        executeTestParallel("testVEGenomicallyAnnotated", spec);
//        //executeTest("testVEGenomicallyAnnotated", spec);
//    }
//
//    @Test
//    public void testVEWriteVCF() {
//        String extraArgs = "-L 1:1-10,000,000 -NO_HEADER -family NA19238+NA19239=NA19240 -mvq 30 -EV MendelianViolationEvaluator";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs + " -o %s -outputVCF %s -NO_HEADER",
//                    2,
//                    Arrays.asList("50321436a65ef7d574286cb0a1c55f7e", "d4bdd06ed5cb1aff1dfee8b69d5d17b8"));
//            executeTestParallel("testVEWriteVCF", spec);
//            //executeTest("testVEWriteVCF", spec);
//        }
//    }

    @Test
    public void testCompVsEvalAC() {
        String extraArgs = "-T VariantEval -R "+b36KGReference+" -o %s -EV GenotypeConcordance -B:evalYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.ug.very.few.lines.vcf -B:compYRI,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk.fake.genotypes.ac.test.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("929e4ec46fb6957c29803531322bb35e"));
        //executeTestParallel("testACDiscordanceAtAC1EvalAC2Comp",spec);
        executeTest("testCompVsEvalAC",spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" -B:eval,vcf " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -tf " + testDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("68044a69f03ba4cc11d2061cc96e9eb5"));
        //executeTestParallel("testTranches",spec);
        executeTest("testTranches",spec);
    }

//    @Test
//    public void testVEValidatePass() {
//        String extraArgs = "-L 1:1-10,000,000";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(withValidateTiTv(withSelect(tests, "DP < 50", "DP50"), 1.0, 4.0) + " " + extraArgs + " -o %s",
//                    1, Arrays.asList("8a0203f0533b628ad7f1f230a43f105f"));
//            executeTestParallel("testVEValidatePass", spec);
//        }
//    }
//
//    @Test
//    public void testVEValidateFail() {
//        String extraArgs = "-L 1:1-10,000,000";
//        for (String tests : testsEnumerations) {
//            WalkerTestSpec spec = new WalkerTestSpec(withValidateTiTv(withSelect(tests, "DP < 50", "DP50"), 1.0, 1.2) + " " + extraArgs + " -o %s",
//                    1, UserException.class);
//            executeTestParallel("testVEValidateFail", spec);
//        }
//    }
//
//    private static String withValidateTiTv(String cmd, double min, double max) {
//        return String.format("%s -validate 'eval.comp_genotypes.all.called.all.titv.tiTvRatio >= %2$s' -validate 'eval.comp_genotypes.all.called.all.titv.tiTvRatio <= %3$s'", cmd, min, max);
//    }
}
