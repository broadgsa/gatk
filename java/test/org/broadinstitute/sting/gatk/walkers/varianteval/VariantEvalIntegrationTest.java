package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta -reportType Grep";

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B eval,VCF," + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    @Test
    public void testVESimple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "bb96e002225df21a84ba7c72613cc67a");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "7427803b4121665e130c76191f862231");

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

        String matchingMD5 = "920940d4006d5d140183f3ba35cafe00";
        expectations.put("", matchingMD5);
        expectations.put(" -known comp_hapmap -known dbsnp", matchingMD5);
        expectations.put(" -known comp_hapmap", "ae9ab5a1556e773b3b1ba80453b14eda");

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
                Arrays.asList("40fca380c2e768cfff5febefc4a73aa0", "a3ce1d70d8ae3874807e9d61994d42af"));
        executeTest("testVEWriteVCF", spec);
    }
}
