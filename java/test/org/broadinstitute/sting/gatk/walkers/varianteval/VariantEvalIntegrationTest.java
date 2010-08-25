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
                    1, Arrays.asList("482c868400f59e17dbc59d667b4b2eca"));
            executeTest("testSelect1", spec);
        }
    }

    @Test
    public void testSelect2() {
        String extraArgs = "-L 1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
                1, Arrays.asList("4dfaa72b23ce297a2c29d9f7e9661c37"));
        executeTest("testSelect2", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFiles[] = {"GenotypeConcordanceEval.vcf", "GenotypeConcordanceEval.vcf.gz"};
        for (String vcfFile : vcfFiles) {
            WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -B:eval,VCF " + validationDataLocation + vcfFile + " -B:comp,VCF " + validationDataLocation + "GenotypeConcordanceComp.vcf -noStandard -E GenotypeConcordance -reportType CSV -o %s",
                    1,
                    Arrays.asList("15d1075d384da2bb7445f7493f2b6a07"));
            executeTest("testVEGenotypeConcordance" + vcfFile, spec);
        }

    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "8891969e7522e728b64c112a2b2f9d1e");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0 -E MendelianViolationEvaluator", "ace2f6170e740a9ee6abc25f130c6848");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();
            for (String tests : testsEnumerations) {
                WalkerTestSpec spec = new WalkerTestSpec( tests + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
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


        String matchingMD5 = "dd513bc72860133a58e9ee542782162b";
        expectations.put("", matchingMD5);
        expectations.put(" -known comp_hapmap -known dbsnp", matchingMD5);
        expectations.put(" -known comp_hapmap", "bef6d1e5fa3a79faf745711e0d8fa2dd");
        for (String tests : testsEnumerations) {
            for (Map.Entry<String, String> entry : expectations.entrySet()) {
                String extraArgs2 = entry.getKey();
                String md5 = entry.getValue();

                WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs1 + extraArgs2 + " -o %s",
                        1, // just one output file
                        Arrays.asList(md5));
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
        executeTest("testVEGenomicallyAnnotated", spec);
    }

    @Test
    public void testVEWriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 30 -E MendelianViolationEvaluator";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs + " -o %s -outputVCF %s",
                    2,
                    Arrays.asList("77abdb58b3166d87daadf397e7fb51c4", "989bc30dea6c8a4cf771cd1b9fdab488"));
            executeTest("testVEWriteVCF", spec);
        }
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }
}
