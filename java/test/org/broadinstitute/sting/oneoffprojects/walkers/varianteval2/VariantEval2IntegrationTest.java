package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class VariantEval2IntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T VariantEval2" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta";

    private static String root = cmdRoot +
            " -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
            " -B eval,VCF,/humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10,000,000", "e7dba09b1856b9be86816939596a5062");
        expectations.put("-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 0", "c170c672ca2ef86068cc5dee9aaac022");
    }

    @Test
    public void testVE2Simple() {

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
    public void testVE2WriteVCF() {
        String extraArgs = "-L 1:1-10,000,000 -family NA19238+NA19239=NA19240 -MVQ 30";
        WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s -outputVCF %s",
                2,
                Arrays.asList("c53d7638df2d7440dee1fd274d1f6384", "9ec81f7389c0971e44e4b8d2d4af3008"));
        executeTest("testVE2WriteVCF", spec);
    }
}