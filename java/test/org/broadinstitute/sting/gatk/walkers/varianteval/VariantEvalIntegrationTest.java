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
        expectations.put("-L 1:1-10,000,000", "5ea9dd109c45c47243f40af759b36d09");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "8b126d445eaec0263d2af6fbfea5ad8b");

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

        expectations.put("", "1bc8126366117045bb08872e098502af");
        expectations.put(" -known comp_hapmap -known dbsnp", "1bc8126366117045bb08872e098502af");
        expectations.put(" -known comp_hapmap", "11e2a3ae77b3ce4b659baa832388f920");

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
                Arrays.asList("cd4813cc4bcce0f8f81d42bcde11f7ce", "a3ce1d70d8ae3874807e9d61994d42af"));
        executeTest("testVEWriteVCF", spec);
    }
}
