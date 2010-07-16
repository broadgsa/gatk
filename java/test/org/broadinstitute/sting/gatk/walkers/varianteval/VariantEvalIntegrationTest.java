package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class
        VariantEvalIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta";

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B eval,VCF," + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -B comp_genotypes,VCF," + validationDataLocation + "yri.trio.gatk.ug.head.vcf -reportType Grep";

    private static String rootGZ = cmdRoot +
                " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
                " -B eval,VCF," + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf.gz" +
                " -B comp_genotypes,VCF," + validationDataLocation + "yri.trio.gatk.ug.head.vcf.gz -reportType Grep";

    private static String[] testsEnumerations = {root, rootGZ};

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -o %s",
                    1, Arrays.asList("119601d7e9e67a1053663b2e546250ed"));
            executeTest("testSelect1", spec);
        }
    }

    @Test
    public void testSelect2() {
        String extraArgs = "-L 1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
                1, Arrays.asList("06d495ab8169a2570eebdc54ecdffe10"));
        executeTest("testSelect2", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        String vcfFiles[] = {"GenotypeConcordanceEval.vcf", "GenotypeConcordanceEval.vcf.gz"};
        for (String vcfFile : vcfFiles) {
            WalkerTestSpec spec = new WalkerTestSpec(cmdRoot + " -B eval,VCF," + validationDataLocation + vcfFile + " -B comp,VCF," + validationDataLocation + "GenotypeConcordanceComp.vcf -E GenotypeConcordance -reportType CSV -o %s",
                    1,
                    Arrays.asList("15d1075d384da2bb7445f7493f2b6a07"));
            executeTest("testVEGenotypeConcordance" + vcfFile, spec);
        }

    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "629b8b124306435ff56b66357354dfbc");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "f51c299d500b347d098c7ab25f54a436");

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
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -MVQ 30" +
                " -B dbsnp_130,dbSNP," + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B comp_hapmap,VCF," + validationDataLocation + "CEU_hapmap_nogt_23.vcf";


        String matchingMD5 = "d01725ce4e46c8fea0855a923c1598fd";
        expectations.put("", matchingMD5);
        expectations.put(" -known comp_hapmap -known dbsnp", matchingMD5);
        expectations.put(" -known comp_hapmap", "a50be9240f6c90503fb6333d8a78b974");
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
        String vecmd = "-T VariantEval " +
                       "-R " + oneKGLocation + "reference/human_b36_both.fasta " +
                       "-L 21 " +
                       "-D " + GATKDataLocation + "dbsnp_129_b36.rod " +
                       "-E CountFunctionalClasses " +
                       "-B eval,VCF," + validationDataLocation + "test.filtered.maf_annotated.vcf " +
                       "-o %s";
        String md5 = "d41d8cd98f00b204e9800998ecf8427e";

        WalkerTestSpec spec = new WalkerTestSpec(vecmd, 1, Arrays.asList(md5));
        executeTest("testVEGenomicallyAnnotated", spec);
    }

    @Test
    public void testVEWriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 30";
        for (String tests : testsEnumerations) {
            WalkerTestSpec spec = new WalkerTestSpec(tests + " " + extraArgs + " -o %s -outputVCF %s",
                    2,
                    Arrays.asList("483f821ce96f4cf571e9bba356c9f325", "d94328f4a5f7c40e95edf2ef13f38ae0"));
            executeTest("testVEWriteVCF", spec);
        }
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }
}
