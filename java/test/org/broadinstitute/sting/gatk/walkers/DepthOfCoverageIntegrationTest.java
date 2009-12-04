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
        expectations.put("-minMAPQ 1", "8b73fad5cce4620907d5da2a985219d5");
        expectations.put("-minMAPQ 100", "1a959892d8ad0523dac2fb097eacb3c2");
        expectations.put("-minDepth 8", "6d549e5a5c4c55420d68e0221a955a0e");
        expectations.put("-minDepth 10", "a367c894e6a48ebb107d2fe004cdfee7");
        expectations.put("-bySample", "93358437153b4d65bdff747e33de1d63");
        expectations.put("-byRG", "777e8427eb4bdad300b23800cb7b0592");
        expectations.put("-histogram", "96f15e1d9d598d48191e20ee84715d46");
        expectations.put("-bases", "baafcb2b90098cad1c5950da9e9932a6");
        expectations.put("-minMAPQ 1 -bySample -byRG -minDepth 8 -histogram -bases", "bf2094b33e0e10fc11a7216bc1097a8b");
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
                "-T DepthOfCoverage -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam -L 1:10,001,890-10,001,895 -o %s",
                1, // just one output file
                Arrays.asList("a332d1539b29dff615b198818a3d4dd1"));
        executeTest("testDepthOfCoverage454", spec);
    }    
}