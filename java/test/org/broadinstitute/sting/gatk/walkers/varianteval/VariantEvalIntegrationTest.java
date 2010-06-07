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


    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec( withSelect(root, "DP < 50", "DP50") + " " + extraArgs + " -o %s",
                1, Arrays.asList("5a330d359b5c7ea0dfa6698b4830db82"));
        executeTest("testSelect1", spec);
    }

    @Test
    public void testSelect2() {
        String extraArgs = "-L 1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec( withSelect(withSelect(root, "DP < 50", "DP50"), "set==\"Intersection\"", "intersection") + " " + extraArgs + " -o %s",
                1, Arrays.asList("e39d6790e4ee8709dfa2eab8598b168e"));
        executeTest("testSelect2", spec);
    }

    @Test
    public void testVEGenotypeConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -B eval,VCF," + validationDataLocation + "GenotypeConcordanceEval.vcf -B comp,VCF," + validationDataLocation + "GenotypeConcordanceComp.vcf -E GenotypeConcordance -reportType CSV -o %s",
                1, // just one output file
                Arrays.asList("608b32189818df7271c546dd5f257033"));
        executeTest("testVEGenotypeConcordance", spec);
    }

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "08b93c2561742dc480534ae16d2d0508");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "65bf329cf1499544d1eb28de786a4561");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testVESimple", spec);
        }
    }

    @Test
    public void testVEComplex() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -MVQ 30" +
                " -B dbsnp_130,dbSNP," + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B comp_hapmap,VCF," + validationDataLocation + "CEU_hapmap_nogt_23.vcf";


        String matchingMD5 = "fc1d4d2aca0667605a6b4d486fcbedf2";
        expectations.put("", matchingMD5);
        expectations.put(" -known comp_hapmap -known dbsnp", matchingMD5);
        expectations.put(" -known comp_hapmap", "fd9cd9970f7e509ca476502869063929");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs2 = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs1 + extraArgs2 + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testVEComplex", spec);
        }
    }

    @Test
    public void testVEWriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 30";
        WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s -outputVCF %s",
                2,
                Arrays.asList("521837758da151b84fca57fd1bb7dad1", "b4a42c90318adc88361691ece50426f2"));
        executeTest("testVEWriteVCF", spec);
    }

    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }
}
