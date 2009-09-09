package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class DepthOfCoverageIntegrationTest extends WalkerTest {
    private static String root = "-L 1:10,164,500-10,164,520 -R /broad/1KG/reference/human_b36_both.fasta -T DepthOfCoverage -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam";
    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-s COMPACT -minMAPQ 1", "a6b2005e37f48545e9604e0f809cf618");
        expectations.put("-s COMPACT", "c3cd44717487ad3c7d6d969583d9bb8c");
        expectations.put("-s COMPACT -ed", "73aacc9febea393a4c3760bd6c70a557");
        expectations.put("-s NONE", "628201adcef01e8bc3e9cdf96b86003b");
        expectations.put("-s NONE -minMAPQ 1", "2425ccbcfdfb9e9a98fcf30c6f2bfcdd");
        expectations.put("-s NONE -minMAPQ 100", "09eaf5a4c811e4ac74de5d179a9ffbe7");
        expectations.put("-s DETAILED -minMAPQ 1", "9059436a57f8c28ef16ee9b2c7cd3ebb");
        expectations.put("-s DETAILED", "fac3bb643f34eb9d90dac0e0358b55ab");
        expectations.put("-s DETAILED -ed", "4d991221aa5e206de11878e88c380de1");
    }

    @Test
    public void testDepthOfCoverage1() {

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testDepthOfCoverage1", spec);
        }
    }

    @Test
    public void testDepthOfCoverage454() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T DepthOfCoverage -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam -L 1:10,001,890-10,001,895 -s DETAILED -o %s",
                1, // just one output file
                Arrays.asList("5c94916ec193ca825c681dd0175c6164"));
        executeTest("testDepthOfCoverage454", spec);
    }    
}