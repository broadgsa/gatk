package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

public class VariantEval2IntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval2" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta -reportType Grep";

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B eval,VCF," + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    @Test
    public void testVE2Simple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "ae10d06b5d30ee227c0e1e18661a18f5");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "a9a682734dbb9a574383a6e7ab385dcc");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testVE2Simple", spec);
        }
    }

    @Test
    public void testVE2Complex() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        String extraArgs1 = "-L " + validationDataLocation + "chr1_b36_pilot3.interval_list -family NA19238+NA19239=NA19240 -MVQ 30" +
                " -B dbsnp_130,dbSNP," + GATKDataLocation + "dbsnp_130_b36.rod" +
                " -B comp_hapmap,VCF," + validationDataLocation + "CEU_hapmap_nogt_23.vcf";

        String eqMD5s = "0fc336da6b77cd7bf1b9c2568a57e94a"; // next two examples should be the same!
        expectations.put("", eqMD5s);
        expectations.put(" -known comp_hapmap -known dbsnp", eqMD5s);
        expectations.put(" -known comp_hapmap", "b6072559ce2d01309639ebb2f6133fa9");

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs2 = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs1 + extraArgs2 + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testVE2Complex", spec);
        }
    }

    @Test
    public void testVE2WriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 30";
        WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s -outputVCF %s",
                2,
                Arrays.asList("b42a8974ebb4354bf172ae5f81922f9e", "a3ce1d70d8ae3874807e9d61994d42af"));
        executeTest("testVE2WriteVCF", spec);
    }
}
