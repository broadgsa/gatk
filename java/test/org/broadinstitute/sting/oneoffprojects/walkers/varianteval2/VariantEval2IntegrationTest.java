package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class VariantEval2IntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval2" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta -reportType Grep";

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B eval,VCF," + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    @Test
    public void testVE2Simple() {
        HashMap<String, String> expectations = new HashMap<String, String>();
        expectations.put("-L 1:1-10,000,000", "44fca2e43dade85ca16fa9e323e4cc9b");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "6910a9a49625a8df2ad797d751be9ed6");

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

        String eqMD5s = "a525597a34cd3fa84c75263dbd026acb"; // next two examples should be the same!
        expectations.put("", eqMD5s);
        expectations.put(" -known comp_hapmap -known dbsnp", eqMD5s);
        expectations.put(" -known comp_hapmap", "35feda6c361408755e7f1dae245462d7");

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
                Arrays.asList("0d2ac5ea23f8cd3237836080c722711b", "a3ce1d70d8ae3874807e9d61994d42af"));
        executeTest("testVE2WriteVCF", spec);
    }
}
